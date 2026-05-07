#!/usr/bin/env python3
"""
generate_pipeline.py — emit an end-to-end Nanopore assembly + annotation
pipeline as a stack of SGE submit scripts wired together with -hold_jid.

Inputs
------
  --pod5-dir      Folder of *.pod5 files (raw signal data)
  --samplesheet   TSV: barcode, sample, genome_length, genus, species,
                  strain, kingdom, reference_proteins
  --output-dir    Where everything (intermediate + final) gets written
  --config        YAML config (see pipeline_config.yaml)

Stages (each emits one SGE script unless noted)
-----------------------------------------------
  01_basecall       dorado basecaller (POD5 -> basecalls.bam)
  02_demux          dorado demux --emit-fastq (per-barcode FASTQs)
  03_rename         flatten demux output to fastq/<sample>.fastq
  04_filter         chopper -q -l on each fastq
  05_subsample      autocycler subsample per sample (4 subsamples)
  06_assemble_*     autocycler helper x 8 assemblers x 4 subsamples,
                    split into N parallel batches (jobs_per_assembly_batch)
  07_organize       group flat *.fasta into per-sample subdirs
  08_compress       autocycler compress per sample
  09_cluster        autocycler cluster per sample
  10_trim_resolve   autocycler trim + resolve on every cluster_*
  11_combine        autocycler combine per sample
  12_collect        gather consensus_assembly.fasta -> consensus/<sample>.fasta
  13_annotate       prokka per sample using samplesheet metadata

Submit everything with:
  bash <output-dir>/submit_all.sh
"""

from __future__ import annotations

import argparse
import csv
import os
import stat
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

try:
    import yaml
except ImportError:
    sys.stderr.write(
        "PyYAML is required. Install with `pip install pyyaml` or "
        "`mamba install pyyaml`.\n"
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# Samplesheet
# ---------------------------------------------------------------------------

REQUIRED_COLS = [
    "barcode",
    "sample",
    "genome_length",
    "genus",
    "species",
    "strain",
    "kingdom",
    "reference_proteins",
]


@dataclass
class Sample:
    barcode: str
    sample: str
    genome_length: int
    genus: str
    species: str
    strain: str
    kingdom: str
    reference_proteins: str  # may be empty -> prokka without --proteins


def parse_samplesheet(path: Path) -> List[Sample]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in REQUIRED_COLS if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(
                f"Samplesheet {path} missing required columns: {missing}\n"
                f"Required: {REQUIRED_COLS}"
            )
        samples: List[Sample] = []
        for row_num, row in enumerate(reader, start=2):
            # Skip blank lines
            if not any((row.get(c) or "").strip() for c in REQUIRED_COLS):
                continue
            try:
                gl = int((row["genome_length"] or "").replace(",", "").strip())
            except ValueError:
                raise SystemExit(
                    f"Row {row_num}: genome_length is not an integer "
                    f"({row['genome_length']!r})"
                )
            samples.append(
                Sample(
                    barcode=row["barcode"].strip(),
                    sample=row["sample"].strip(),
                    genome_length=gl,
                    genus=row["genus"].strip(),
                    species=row["species"].strip(),
                    strain=row["strain"].strip(),
                    kingdom=row["kingdom"].strip(),
                    reference_proteins=(row.get("reference_proteins") or "").strip(),
                )
            )
    if not samples:
        raise SystemExit(f"Samplesheet {path} has no data rows.")

    # Sanity: unique barcodes and unique sample names
    barcodes = [s.barcode for s in samples]
    if len(set(barcodes)) != len(barcodes):
        raise SystemExit(f"Duplicate barcodes in samplesheet: {barcodes}")
    sample_names = [s.sample for s in samples]
    if len(set(sample_names)) != len(sample_names):
        raise SystemExit(f"Duplicate sample names in samplesheet: {sample_names}")
    return samples


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


def load_config(path: Path) -> dict:
    with path.open() as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise SystemExit(f"Config {path} did not parse to a mapping.")
    return cfg


# ---------------------------------------------------------------------------
# SGE script header
# ---------------------------------------------------------------------------


def sge_header(
    job_name: str,
    stdout_dir: Path,
    h_rt: str,
    h_vmem: str,
    pe_parallel: int,
    email: str,
    *,
    gpu: bool = False,
    gpumem: str = "20G",
    hold_jid: str | None = None,
) -> str:
    lines = [
        "#!/bin/bash",
        f"#$ -N {job_name}",
        f"#$ -o {stdout_dir}/{job_name}.log",
        "#$ -j y",
        "#$ -cwd",
        f"#$ -l h_rt={h_rt}",
        f"#$ -l h_vmem={h_vmem}",
        f"#$ -pe parallel {pe_parallel}",
        "#$ -m a",
        f"#$ -M {email}",
    ]
    if gpu:
        lines += [f"#$ -l gpumem={gpumem}", "#$ -l gpu=1"]
    if hold_jid:
        lines.append(f"#$ -hold_jid {hold_jid}")
    return "\n".join(lines) + "\n"


def cuda_block(cfg: dict) -> str:
    return textwrap.dedent(
        f"""
        # CUDA + NVIDIA driver location
        CUDA_DIR={cfg['cuda_dir']}
        export PATH=$PATH:$CUDA_DIR/{cfg['cuda_version']}/bin
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/nvidia:$CUDA_DIR/{cfg['cuda_version']}/lib64:$CUDA_DIR/{cfg['cuda_version']}/extras/CUPTI/lib64:$CUDA_DIR/cuda-drivers
        """
    ).strip()


def conda_block(env: str) -> str:
    return textwrap.dedent(
        f"""
        # Activate conda env
        . ~/.bashrc
        mamba activate {env}
        """
    ).strip()


# ---------------------------------------------------------------------------
# Stage builders
# ---------------------------------------------------------------------------


@dataclass
class Layout:
    """Resolved output paths."""

    root: Path
    pod5_dir: Path
    basecalled: Path
    bam_path: Path
    fastq_raw: Path
    fastq: Path
    fastq_filtered: Path
    autocycler: Path
    assemblies_flat: Path
    assemblies_fasta: Path
    cluster_dir: Path
    consensus: Path
    prokka: Path
    submit_dir: Path
    stdout_dir: Path

    @classmethod
    def from_root(cls, root: Path, pod5_dir: Path) -> "Layout":
        return cls(
            root=root,
            pod5_dir=pod5_dir,
            basecalled=root / "basecalled",
            bam_path=root / "basecalled" / "basecalls.bam",
            fastq_raw=root / "fastq_raw",
            fastq=root / "fastq",
            fastq_filtered=root / "fastq_filtered",
            autocycler=root / "autocycler",
            assemblies_flat=root / "autocycler" / "assemblies",
            assemblies_fasta=root / "autocycler" / "assemblies_fasta",
            cluster_dir=root / "autocycler" / "cluster",
            consensus=root / "consensus",
            prokka=root / "prokka",
            submit_dir=root / "submit_scripts",
            stdout_dir=root / "stdout",
        )


def write_script(path: Path, body: str) -> None:
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


# ---- 01 basecall ---------------------------------------------------------

def stage_basecall(cfg: dict, layout: Layout) -> str:
    name = "s01_basecall"
    mod_models = ",".join(cfg["mod_models"])
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_gpu"],
        cfg["sge_h_vmem_gpu"],
        cfg["sge_pe_parallel"],
        cfg["email"],
        gpu=True,
        gpumem=cfg["sge_gpumem"],
    )
    body += conda_block(cfg["conda_env_basecall"]) + "\n"
    body += cuda_block(cfg) + "\n\n"
    body += textwrap.dedent(
        f"""
        mkdir -p {layout.basecalled}
        {cfg['dorado_binary']} basecaller \\
            {cfg['basecall_model']} \\
            {layout.pod5_dir} \\
            --modified-bases-models {mod_models} \\
            --kit-name {cfg['kit_name']} \\
            > {layout.bam_path}
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 02 demux ------------------------------------------------------------

def stage_demux(cfg: dict, layout: Layout, hold: str) -> str:
    name = "s02_demux"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_gpu"],
        cfg["sge_h_vmem_gpu"],
        cfg["sge_pe_parallel"],
        cfg["email"],
        gpu=True,
        gpumem=cfg["sge_gpumem"],
        hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_demux"]) + "\n"
    body += cuda_block(cfg) + "\n\n"
    body += textwrap.dedent(
        f"""
        mkdir -p {layout.fastq_raw}
        {cfg['dorado_binary']} demux \\
            --output-dir {layout.fastq_raw} \\
            --no-classify \\
            --emit-fastq \\
            --threads {cfg['sge_pe_parallel']} \\
            --emit-summary \\
            {layout.bam_path}
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 03 rename -----------------------------------------------------------

def stage_rename(
    cfg: dict, layout: Layout, samples: List[Sample], hold: str
) -> str:
    """Flatten dorado demux output into fastq/<sample>.fastq.

    dorado demux writes files containing the barcode in the filename
    (e.g. <prefix>_barcode01.fastq). We rglob for files matching each
    barcode and concatenate into <sample>.fastq.
    """
    name = "s03_rename"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_cpu"],
        cfg["sge_h_vmem_cpu"],
        1,
        cfg["email"],
        hold_jid=hold,
    )
    body += textwrap.dedent(
        f"""

        set -euo pipefail
        mkdir -p {layout.fastq}

        """
    )
    for s in samples:
        body += textwrap.dedent(
            f"""
            # {s.sample} <- {s.barcode}
            files=$(find {layout.fastq_raw} -type f \\( -name "*{s.barcode}*.fastq" -o -name "*{s.barcode}*.fastq.gz" \\))
            if [ -z "$files" ]; then
                echo "ERROR: no demux output found for barcode {s.barcode} (sample {s.sample})" >&2
                exit 1
            fi
            : > {layout.fastq}/{s.sample}.fastq
            for f in $files; do
                case "$f" in
                    *.gz) zcat "$f" >> {layout.fastq}/{s.sample}.fastq ;;
                    *)    cat  "$f" >> {layout.fastq}/{s.sample}.fastq ;;
                esac
            done
            echo "Wrote {layout.fastq}/{s.sample}.fastq ($(wc -l < {layout.fastq}/{s.sample}.fastq) lines)"

            """
        )
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 04 filter (chopper) -------------------------------------------------

def stage_filter(
    cfg: dict, layout: Layout, samples: List[Sample], hold: str
) -> str:
    name = "s04_filter"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_cpu"],
        cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"],
        cfg["email"],
        hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_chopper"]) + "\n\n"
    body += textwrap.dedent(
        f"""
        set -euo pipefail
        mkdir -p {layout.fastq_filtered}

        """
    )
    for s in samples:
        body += textwrap.dedent(
            f"""
            echo "chopper {s.sample}..."
            chopper -q {cfg['chopper_quality']} -l {cfg['chopper_min_length']} -t {cfg['sge_pe_parallel']} \\
                -i {layout.fastq}/{s.sample}.fastq \\
                > {layout.fastq_filtered}/{s.sample}.fastq

            """
        )
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 05 autocycler subsample --------------------------------------------

def stage_subsample(
    cfg: dict, layout: Layout, samples: List[Sample], hold: str
) -> str:
    name = "s05_subsample"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_cpu"],
        cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"],
        cfg["email"],
        hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
    body += f"set -euo pipefail\nmkdir -p {layout.autocycler}\n\n"
    for s in samples:
        body += textwrap.dedent(
            f"""
            echo "autocycler subsample {s.sample} (genome_size={s.genome_length})"
            autocycler subsample \\
                --reads {layout.fastq_filtered}/{s.sample}.fastq \\
                --out_dir {layout.autocycler}/{s.sample} \\
                --genome_size {s.genome_length}

            """
        )
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 06 autocycler helper (assemblies) ----------------------------------

def stage_assemble(
    cfg: dict, layout: Layout, samples: List[Sample], hold: str
) -> List[str]:
    """Build all helper commands across samples x assemblers x subsamples,
    then split into batches of `jobs_per_assembly_batch` per SGE script.
    """
    assemblers: List[str] = cfg["assemblers"]
    threads_default = cfg["threads_per_assembler"]
    threads_necat = cfg["threads_necat"]
    jobs_per = cfg["jobs_per_assembly_batch"]
    # Subsample numbers used by autocycler subsample
    subsample_nums = ["01", "02", "03", "04"]

    commands: List[str] = []
    for s in samples:
        for sub in subsample_nums:
            reads = f"{layout.autocycler}/{s.sample}/sample_{sub}.fastq"
            for asm in assemblers:
                t = threads_necat if asm == "necat" else threads_default
                out_prefix = (
                    f"{layout.assemblies_flat}/{s.sample}_{asm}_{sub}"
                )
                commands.append(
                    f"autocycler helper {asm} "
                    f"--reads {reads} "
                    f"--out_prefix {out_prefix} "
                    f"--threads {t} "
                    f"--genome_size {s.genome_length}"
                )

    # Chunk into jobs_per-sized batches
    batches = [commands[i : i + jobs_per] for i in range(0, len(commands), jobs_per)]
    names: List[str] = []
    for i, batch in enumerate(batches, start=1):
        name = f"s06_assemble_batch_{i:03d}"
        body = sge_header(
            name,
            layout.stdout_dir,
            cfg["sge_h_rt_assembly"],
            cfg["sge_h_vmem_cpu"],
            cfg["sge_pe_parallel"],
            cfg["email"],
            hold_jid=hold,
        )
        body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
        body += "set -euo pipefail\n"
        body += f"mkdir -p {layout.assemblies_flat}\n\n"
        body += "\n".join(batch) + "\n"
        write_script(layout.submit_dir / f"{name}.sh", body)
        names.append(name)
    return names


# ---- 07 organize ---------------------------------------------------------

def stage_organize(cfg: dict, layout: Layout, hold: str) -> str:
    """Group flat assembly fastas into per-sample subdirs by the first
    three underscore-delimited tokens of the filename.
    """
    name = "s07_organize"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_cpu"],
        cfg["sge_h_vmem_cpu"],
        1,
        cfg["email"],
        hold_jid=hold,
    )
    body += textwrap.dedent(
        f"""

        set -euo pipefail
        python3 - <<'PYEOF'
        import shutil
        from collections import defaultdict
        from pathlib import Path

        src = Path("{layout.assemblies_flat}")
        dst = Path("{layout.assemblies_fasta}")
        dst.mkdir(parents=True, exist_ok=True)

        files = list(src.glob("*.fasta"))
        groups = defaultdict(list)
        for f in files:
            parts = f.name.split("_")
            # autocycler helper out_prefix is <sample>_<asm>_<sub>; the
            # produced fasta filename includes one or more dot-suffixes.
            # We group by everything up to the assembler token.
            # Find index of first known assembler token.
            assemblers = {{"canu", "flye", "metamdbg", "miniasm", "necat",
                          "nextdenovo", "plassembler", "raven"}}
            idx = None
            for i, p in enumerate(parts):
                if p in assemblers:
                    idx = i
                    break
            if idx is None:
                # fall back to first three tokens
                idx = 3
            sample = "_".join(parts[:idx]) if idx else f.stem
            groups[sample].append(f)

        for sample, members in groups.items():
            sub = dst / sample
            sub.mkdir(exist_ok=True)
            for f in members:
                shutil.copy2(f, sub / f.name)
            print(f"{{sample}}: {{len(members)}} assemblies")
        PYEOF
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 08 compress ---------------------------------------------------------

def stage_compress(
    cfg: dict, layout: Layout, samples: List[Sample], hold: str
) -> str:
    name = "s08_compress"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_cpu"],
        cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"],
        cfg["email"],
        hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
    body += "set -euo pipefail\n\n"
    for s in samples:
        body += textwrap.dedent(
            f"""
            echo "autocycler compress {s.sample}"
            autocycler compress \\
                --threads {cfg['sge_pe_parallel']} \\
                --assemblies_dir {layout.assemblies_fasta}/{s.sample} \\
                --autocycler_dir {layout.autocycler}/{s.sample}

            """
        )
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 09 cluster ----------------------------------------------------------

def stage_cluster(
    cfg: dict, layout: Layout, samples: List[Sample], hold: str
) -> str:
    name = "s09_cluster"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_cpu"],
        cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"],
        cfg["email"],
        hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
    body += "set -euo pipefail\n\n"
    body += f"mkdir -p {layout.cluster_dir}\n\n"
    for s in samples:
        # The original moves the autocycler output for a sample under cluster/<sample>
        # before running cluster. We mirror that.
        body += textwrap.dedent(
            f"""
            echo "autocycler cluster {s.sample}"
            mkdir -p {layout.cluster_dir}/{s.sample}
            # cluster works in-place on the autocycler dir
            autocycler cluster \\
                --autocycler_dir {layout.autocycler}/{s.sample}
            # mirror the layout the trim/resolve/combine steps expect
            ln -sfn {layout.autocycler}/{s.sample}/clustering {layout.cluster_dir}/{s.sample}/clustering

            """
        )
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 10 trim + resolve --------------------------------------------------

def stage_trim_resolve(cfg: dict, layout: Layout, hold: str) -> str:
    name = "s10_trim_resolve"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_cpu"],
        cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"],
        cfg["email"],
        hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
    body += textwrap.dedent(
        f"""
        set -euo pipefail

        find {layout.cluster_dir} -type d -path "*/clustering/qc_pass" | while read qc_dir; do
            for c in "$qc_dir"/cluster_*; do
                if [ -d "$c" ]; then
                    echo "trim $c"
                    autocycler trim -c "$c"
                fi
            done
        done

        find {layout.cluster_dir} -type d -path "*/clustering/qc_pass" | while read qc_dir; do
            for c in "$qc_dir"/cluster_*; do
                if [ -d "$c" ]; then
                    echo "resolve $c"
                    autocycler resolve -c "$c"
                fi
            done
        done
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 11 combine ---------------------------------------------------------

def stage_combine(
    cfg: dict, layout: Layout, samples: List[Sample], hold: str
) -> str:
    name = "s11_combine"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_cpu"],
        cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"],
        cfg["email"],
        hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
    body += "set -euo pipefail\n\n"
    for s in samples:
        body += textwrap.dedent(
            f"""
            sample_dir={layout.cluster_dir}/{s.sample}
            gfas=$(ls $sample_dir/clustering/qc_pass/cluster_*/5_final.gfa 2>/dev/null || true)
            if [ -z "$gfas" ]; then
                echo "WARNING: no 5_final.gfa for {s.sample}, skipping combine"
                continue
            fi
            echo "autocycler combine {s.sample}"
            autocycler combine -a $sample_dir -i $gfas

            """
        )
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 12 collect ---------------------------------------------------------

def stage_collect(
    cfg: dict, layout: Layout, samples: List[Sample], hold: str
) -> str:
    name = "s12_collect"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_cpu"],
        cfg["sge_h_vmem_cpu"],
        1,
        cfg["email"],
        hold_jid=hold,
    )
    body += textwrap.dedent(
        f"""

        set -euo pipefail
        mkdir -p {layout.consensus}

        """
    )
    for s in samples:
        body += textwrap.dedent(
            f"""
            src={layout.cluster_dir}/{s.sample}/consensus_assembly.fasta
            if [ ! -f "$src" ]; then
                # fall back to autocycler dir if combine wrote there instead
                src={layout.autocycler}/{s.sample}/consensus_assembly.fasta
            fi
            if [ -f "$src" ]; then
                cp "$src" {layout.consensus}/{s.sample}.fasta
                echo "Collected {s.sample} <- $src"
            else
                echo "ERROR: no consensus_assembly.fasta found for {s.sample}" >&2
                exit 1
            fi

            """
        )
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- 13 prokka ----------------------------------------------------------

def stage_annotate(
    cfg: dict, layout: Layout, samples: List[Sample], hold: str
) -> str:
    name = "s13_annotate"
    body = sge_header(
        name,
        layout.stdout_dir,
        cfg["sge_h_rt_cpu"],
        cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"],
        cfg["email"],
        hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_prokka"]) + "\n\n"
    body += f"set -euo pipefail\nmkdir -p {layout.prokka}\n\n"
    for s in samples:
        proteins = f"--proteins {s.reference_proteins}" if s.reference_proteins else ""
        body += textwrap.dedent(
            f"""
            echo "prokka {s.sample}"
            prokka --dbdir {cfg['prokka_db_dir']} \\
                --outdir {layout.prokka}/{s.sample} \\
                --prefix {s.sample} \\
                --locustag {s.sample} \\
                --genus {s.genus} \\
                --species {s.species} \\
                --strain {s.strain} \\
                --kingdom {s.kingdom} \\
                {proteins} \\
                --cdsrnaolap \\
                --cpus {cfg['prokka_threads']} \\
                --rfam \\
                {layout.consensus}/{s.sample}.fasta

            """
        )
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def write_submit_all(
    layout: Layout, ordered_jobs: List[List[str]]
) -> None:
    """ordered_jobs is a list of stages; each stage is a list of job names
    that should be submitted in parallel and held by the next stage."""
    lines = [
        "#!/bin/bash",
        "# Auto-generated by generate_pipeline.py",
        "# Submits every stage with -hold_jid wired to the previous stage.",
        "set -euo pipefail",
        f"cd {layout.submit_dir}",
        "",
    ]
    for stage in ordered_jobs:
        for job in stage:
            lines.append(f"qsub {job}.sh")
        lines.append("")
    body = "\n".join(lines) + "\n"
    out = layout.root / "submit_all.sh"
    out.write_text(body)
    out.chmod(out.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pod5-dir", required=True, type=Path, help="Folder of *.pod5 files")
    p.add_argument("--samplesheet", required=True, type=Path, help="TSV samplesheet")
    p.add_argument("--output-dir", required=True, type=Path, help="Where outputs go")
    p.add_argument("--config", required=True, type=Path, help="YAML config")
    args = p.parse_args()

    if not args.pod5_dir.is_dir():
        raise SystemExit(f"--pod5-dir does not exist or is not a directory: {args.pod5_dir}")
    if not args.samplesheet.is_file():
        raise SystemExit(f"--samplesheet not found: {args.samplesheet}")
    if not args.config.is_file():
        raise SystemExit(f"--config not found: {args.config}")

    cfg = load_config(args.config)
    samples = parse_samplesheet(args.samplesheet)
    layout = Layout.from_root(args.output_dir.resolve(), args.pod5_dir.resolve())

    # Make all the output subdirs (so SGE -o paths exist before submission)
    for d in [
        layout.root,
        layout.submit_dir,
        layout.stdout_dir,
        layout.basecalled,
        layout.fastq_raw,
        layout.fastq,
        layout.fastq_filtered,
        layout.autocycler,
        layout.assemblies_flat,
        layout.assemblies_fasta,
        layout.cluster_dir,
        layout.consensus,
        layout.prokka,
    ]:
        d.mkdir(parents=True, exist_ok=True)

    # Build stages
    n01 = stage_basecall(cfg, layout)
    n02 = stage_demux(cfg, layout, hold=n01)
    n03 = stage_rename(cfg, layout, samples, hold=n02)
    n04 = stage_filter(cfg, layout, samples, hold=n03)
    n05 = stage_subsample(cfg, layout, samples, hold=n04)
    n06_list = stage_assemble(cfg, layout, samples, hold=n05)
    # The next stage needs to wait on every batch -> SGE accepts comma-separated
    n06_hold = ",".join(n06_list)
    n07 = stage_organize(cfg, layout, hold=n06_hold)
    n08 = stage_compress(cfg, layout, samples, hold=n07)
    n09 = stage_cluster(cfg, layout, samples, hold=n08)
    n10 = stage_trim_resolve(cfg, layout, hold=n09)
    n11 = stage_combine(cfg, layout, samples, hold=n10)
    n12 = stage_collect(cfg, layout, samples, hold=n11)
    n13 = stage_annotate(cfg, layout, samples, hold=n12)

    ordered = [
        [n01], [n02], [n03], [n04], [n05],
        n06_list,
        [n07], [n08], [n09], [n10], [n11], [n12], [n13],
    ]
    write_submit_all(layout, ordered)

    # Summary
    print(f"Generated pipeline at: {layout.root}")
    print(f"  - submit_scripts/   ({sum(len(s) for s in ordered)} scripts)")
    print(f"  - submit_all.sh     (driver)")
    print(f"  - {len(samples)} sample(s): {[s.sample for s in samples]}")
    print()
    print(f"Run with:")
    print(f"  bash {layout.root}/submit_all.sh")


if __name__ == "__main__":
    main()
