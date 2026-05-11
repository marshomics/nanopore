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
  s01_basecall       dorado basecaller (POD5 -> basecalls.bam)
  s02_demux          dorado demux --emit-fastq (per-barcode FASTQs)
  s03_rename         flatten demux output to fastq/<sample>.fastq
  s04_filter         chopper -q -l on each fastq
  s05_subsample      autocycler subsample per sample (4 subsamples)
  s06_assemble_*     autocycler helper x 8 assemblers x 4 subsamples,
                     split into N parallel batches (jobs_per_assembly_batch)
  s07_organize       group flat *.fasta into per-sample subdirs AND drop
                     contigs below min_contig_length (fixes compress's
                     "contigs per input exceeds threshold" error)
  s08_compress       autocycler compress per sample
  s09_cluster        autocycler cluster per sample
  s10_trim_resolve   autocycler trim + resolve on every cluster_*
  s11_combine        autocycler combine per sample
  s12_collect        gather consensus_assembly.fasta -> consensus/<sample>.fasta
  s13_annotate       prokka per sample using samplesheet metadata

Logging + resume
----------------
Every stage appends START / SKIP / OK / WARN / FAIL / DONE events to
<output-dir>/pipeline.log so you can `tail -f` a single file across all
stages. Each per-sample step also checks whether its expected output
already exists and skips if so; re-running submit_all.sh after a partial
failure picks up where the previous run stopped.

To force a stage to re-run a sample, delete that sample's output (e.g.
rm -r <output>/autocycler/<sample>) and resubmit.

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
from dataclasses import dataclass
from pathlib import Path
from typing import List

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
    cfg.setdefault("min_contig_length", 1000)
    cfg.setdefault("max_contigs", 25)
    return cfg


# ---------------------------------------------------------------------------
# SGE script header + shared bash helpers
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
        CUDA_DIR={cfg['cuda_dir']}
        export PATH=$PATH:$CUDA_DIR/{cfg['cuda_version']}/bin
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/nvidia:$CUDA_DIR/{cfg['cuda_version']}/lib64:$CUDA_DIR/{cfg['cuda_version']}/extras/CUPTI/lib64:$CUDA_DIR/cuda-drivers
        """
    ).strip()


def conda_block(env: str) -> str:
    return textwrap.dedent(
        f"""
        . ~/.bashrc
        mamba activate {env}
        """
    ).strip()


def log_setup(name: str, layout: "Layout") -> str:
    """Bash block injected near the top of every stage script. Wires up
    the shared pipeline.log, START/DONE/FAIL trap, and a `log_event`
    helper that scripts call for SKIP/OK/WARN per-sample events."""
    return textwrap.dedent(
        f"""
        PIPELINE_LOG={layout.root}/pipeline.log
        STAGE_LOG={layout.stdout_dir}/{name}.log
        JOB_NAME={name}
        STAGE_FAIL_COUNT=0
        log_event() {{
            local level="$1"; shift
            local ts
            ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
            printf '%s [%s] %s %s\\n' "$ts" "$level" "$JOB_NAME" "$*" >> "$PIPELINE_LOG"
            printf '%s [%s] %s %s\\n' "$ts" "$level" "$JOB_NAME" "$*" >&2
        }}
        log_warn_errors_from() {{
            local f="$1"
            [ -f "$f" ] || return 0
            grep -niE 'error|warn|fail|traceback' "$f" 2>/dev/null \\
                | head -50 \\
                | while IFS= read -r line; do
                    log_event SCAN "$line"
                done
        }}
        on_exit() {{
            local rc=$?
            if [ "$rc" -eq 0 ] && [ "$STAGE_FAIL_COUNT" -eq 0 ]; then
                log_event DONE "stage finished successfully"
            elif [ "$rc" -eq 0 ]; then
                log_event DONE "stage finished with $STAGE_FAIL_COUNT per-sample failures"
            else
                log_event FAIL "stage exited with code $rc"
            fi
        }}
        trap on_exit EXIT
        log_event START
        """
    ).strip()


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


@dataclass
class Layout:
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


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------


# ---- s01 basecall --------------------------------------------------------

def stage_basecall(cfg: dict, layout: Layout) -> str:
    name = "s01_basecall"
    mod_models = ",".join(cfg["mod_models"])
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_gpu"], cfg["sge_h_vmem_gpu"],
        cfg["sge_pe_parallel"], cfg["email"], gpu=True, gpumem=cfg["sge_gpumem"],
    )
    body += conda_block(cfg["conda_env_basecall"]) + "\n"
    body += cuda_block(cfg) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += textwrap.dedent(
        f"""
        mkdir -p {layout.basecalled}

        if [ -s {layout.bam_path} ]; then
            log_event SKIP "basecalls.bam already exists ($(stat -c %s {layout.bam_path}) bytes)"
        else
            log_event OK "starting dorado basecaller"
            set +e
            {cfg['dorado_binary']} basecaller \\
                {cfg['basecall_model']} \\
                {layout.pod5_dir} \\
                --modified-bases-models {mod_models} \\
                --kit-name {cfg['kit_name']} \\
                > {layout.bam_path}
            rc=$?
            set -e
            if [ $rc -ne 0 ]; then
                log_event WARN "dorado basecaller exited with code $rc"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                rm -f {layout.bam_path}
                exit $rc
            fi
            log_event OK "basecalls.bam written"
        fi

        log_warn_errors_from "$STAGE_LOG"
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s02 demux -----------------------------------------------------------

def stage_demux(cfg: dict, layout: Layout, hold: str) -> str:
    name = "s02_demux"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_gpu"], cfg["sge_h_vmem_gpu"],
        cfg["sge_pe_parallel"], cfg["email"], gpu=True, gpumem=cfg["sge_gpumem"],
        hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_demux"]) + "\n"
    body += cuda_block(cfg) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += textwrap.dedent(
        f"""
        mkdir -p {layout.fastq_raw}

        # Resume if any fastq output already exists in the demux dir.
        if find {layout.fastq_raw} -type f \\( -name "*.fastq" -o -name "*.fastq.gz" \\) | grep -q .; then
            log_event SKIP "demux output already present in {layout.fastq_raw}"
        else
            log_event OK "starting dorado demux"
            set +e
            {cfg['dorado_binary']} demux \\
                --output-dir {layout.fastq_raw} \\
                --no-classify \\
                --emit-fastq \\
                --threads {cfg['sge_pe_parallel']} \\
                --emit-summary \\
                {layout.bam_path}
            rc=$?
            set -e
            if [ $rc -ne 0 ]; then
                log_event WARN "dorado demux exited with code $rc"
                exit $rc
            fi
            log_event OK "demux complete"
        fi

        log_warn_errors_from "$STAGE_LOG"
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s03 rename ----------------------------------------------------------

def stage_rename(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s03_rename"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        1, cfg["email"], hold_jid=hold,
    )
    body += log_setup(name, layout) + "\n\n"
    body += textwrap.dedent(
        f"""
        mkdir -p {layout.fastq}

        """
    )
    for s in samples:
        out = f"{layout.fastq}/{s.sample}.fastq"
        body += textwrap.dedent(
            f"""
            # {s.sample} <- {s.barcode}
            if [ -s {out} ]; then
                log_event SKIP "{s.sample} (fastq exists)"
            else
                files=$(find {layout.fastq_raw} -type f \\( -name "*{s.barcode}*.fastq" -o -name "*{s.barcode}*.fastq.gz" \\) 2>/dev/null)
                if [ -z "$files" ]; then
                    log_event WARN "{s.sample}: no demux output for {s.barcode}"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    : > {out}
                    for f in $files; do
                        case "$f" in
                            *.gz) zcat "$f" >> {out} ;;
                            *)    cat  "$f" >> {out} ;;
                        esac
                    done
                    nlines=$(wc -l < {out})
                    log_event OK "{s.sample} wrote $nlines lines"
                fi
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s04 filter ----------------------------------------------------------

def stage_filter(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s04_filter"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_chopper"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += f"mkdir -p {layout.fastq_filtered}\n\n"
    for s in samples:
        inp = f"{layout.fastq}/{s.sample}.fastq"
        out = f"{layout.fastq_filtered}/{s.sample}.fastq"
        body += textwrap.dedent(
            f"""
            if [ -s {out} ]; then
                log_event SKIP "{s.sample} (filtered fastq exists)"
            elif [ ! -s {inp} ]; then
                log_event WARN "{s.sample}: input {inp} missing/empty"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                set +e
                chopper -q {cfg['chopper_quality']} -l {cfg['chopper_min_length']} -t {cfg['sge_pe_parallel']} \\
                    -i {inp} > {out}
                rc=$?
                set -e
                if [ $rc -ne 0 ]; then
                    log_event WARN "{s.sample}: chopper rc=$rc"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    rm -f {out}
                else
                    log_event OK "{s.sample} chopper done"
                fi
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s05 subsample -------------------------------------------------------

def stage_subsample(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s05_subsample"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += f"mkdir -p {layout.autocycler}\n\n"
    for s in samples:
        sample_dir = f"{layout.autocycler}/{s.sample}"
        body += textwrap.dedent(
            f"""
            if [ -f {sample_dir}/sample_04.fastq ]; then
                log_event SKIP "{s.sample} (4 subsamples already exist)"
            else
                set +e
                autocycler subsample \\
                    --reads {layout.fastq_filtered}/{s.sample}.fastq \\
                    --out_dir {sample_dir} \\
                    --genome_size {s.genome_length}
                rc=$?
                set -e
                if [ $rc -ne 0 ]; then
                    log_event WARN "{s.sample}: subsample rc=$rc"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    log_event OK "{s.sample} subsampled (genome_size={s.genome_length})"
                fi
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s06 assemble --------------------------------------------------------

def stage_assemble(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> List[str]:
    assemblers: List[str] = cfg["assemblers"]
    threads_default = cfg["threads_per_assembler"]
    threads_necat = cfg["threads_necat"]
    jobs_per = cfg["jobs_per_assembly_batch"]
    subsample_nums = ["01", "02", "03", "04"]

    commands: List[tuple] = []  # (out_prefix, cmd_line, label)
    for s in samples:
        for sub in subsample_nums:
            reads = f"{layout.autocycler}/{s.sample}/sample_{sub}.fastq"
            for asm in assemblers:
                t = threads_necat if asm == "necat" else threads_default
                out_prefix = f"{layout.assemblies_flat}/{s.sample}_{asm}_{sub}"
                cmd = (
                    f"autocycler helper {asm} "
                    f"--reads {reads} "
                    f"--out_prefix {out_prefix} "
                    f"--threads {t} "
                    f"--genome_size {s.genome_length}"
                )
                label = f"{s.sample}_{asm}_{sub}"
                commands.append((out_prefix, cmd, label))

    batches = [commands[i : i + jobs_per] for i in range(0, len(commands), jobs_per)]
    names: List[str] = []
    for i, batch in enumerate(batches, start=1):
        name = f"s06_assemble_batch_{i:03d}"
        body = sge_header(
            name, layout.stdout_dir, cfg["sge_h_rt_assembly"], cfg["sge_h_vmem_cpu"],
            cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
        )
        body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
        body += log_setup(name, layout) + "\n\n"
        body += f"mkdir -p {layout.assemblies_flat}\n\n"
        for out_prefix, cmd, label in batch:
            body += textwrap.dedent(
                f"""
                if [ -s {out_prefix}.fasta ]; then
                    log_event SKIP "{label} ({out_prefix}.fasta exists)"
                else
                    log_event OK "{label} starting"
                    set +e
                    {cmd}
                    rc=$?
                    set -e
                    if [ $rc -ne 0 ] || [ ! -s {out_prefix}.fasta ]; then
                        log_event WARN "{label} rc=$rc (assembler failures are common, continuing)"
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    else
                        log_event OK "{label} done"
                    fi
                fi

                """
            )
        body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
        write_script(layout.submit_dir / f"{name}.sh", body)
        names.append(name)
    return names


# ---- s07 organize (with contig-length filter) ---------------------------

def stage_organize(cfg: dict, layout: Layout, hold: str) -> str:
    name = "s07_organize"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        1, cfg["email"], hold_jid=hold,
    )
    body += log_setup(name, layout) + "\n\n"
    min_len = int(cfg["min_contig_length"])
    body += textwrap.dedent(
        f"""
        log_event OK "organizing and filtering (min_contig_length={min_len})"
        python3 - <<'PYEOF'
        import os, sys, shutil
        from collections import defaultdict
        from pathlib import Path

        src = Path("{layout.assemblies_flat}")
        dst = Path("{layout.assemblies_fasta}")
        dst.mkdir(parents=True, exist_ok=True)
        MIN_LEN = {min_len}
        PIPELINE_LOG = Path("{layout.root}/pipeline.log")

        def log(level, msg):
            import datetime
            ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
            line = f"{{ts}} [{{level}}] s07_organize {{msg}}\\n"
            with PIPELINE_LOG.open("a") as fh:
                fh.write(line)
            sys.stderr.write(line)

        ASSEMBLERS = {{"canu","flye","metamdbg","miniasm","necat",
                       "nextdenovo","plassembler","raven"}}

        def filter_fasta(in_path: Path, out_path: Path, min_len: int):
            kept = 0
            dropped = 0
            with in_path.open() as fin, out_path.open("w") as fout:
                header = None
                seq = []
                def flush():
                    nonlocal kept, dropped
                    if header is None:
                        return
                    s = "".join(seq)
                    if len(s) >= min_len:
                        fout.write(header + "\\n")
                        for i in range(0, len(s), 80):
                            fout.write(s[i:i+80] + "\\n")
                        kept_local = 1
                        return 1, 0
                    else:
                        return 0, 1
                for line in fin:
                    line = line.rstrip("\\n")
                    if line.startswith(">"):
                        if header is not None:
                            s = "".join(seq)
                            if len(s) >= min_len:
                                fout.write(header + "\\n")
                                for i in range(0, len(s), 80):
                                    fout.write(s[i:i+80] + "\\n")
                                kept += 1
                            else:
                                dropped += 1
                        header = line
                        seq = []
                    else:
                        seq.append(line)
                if header is not None:
                    s = "".join(seq)
                    if len(s) >= min_len:
                        fout.write(header + "\\n")
                        for i in range(0, len(s), 80):
                            fout.write(s[i:i+80] + "\\n")
                        kept += 1
                    else:
                        dropped += 1
            return kept, dropped

        files = list(src.glob("*.fasta"))
        groups = defaultdict(list)
        for f in files:
            parts = f.stem.split("_")
            idx = None
            for i, p in enumerate(parts):
                if p in ASSEMBLERS:
                    idx = i
                    break
            if idx is None:
                idx = 3 if len(parts) >= 3 else len(parts)
            sample = "_".join(parts[:idx]) if idx else f.stem
            groups[sample].append(f)

        total_kept = 0
        total_dropped = 0
        for sample, members in groups.items():
            sub = dst / sample
            sub.mkdir(exist_ok=True)
            for f in members:
                out = sub / f.name
                if out.exists() and out.stat().st_size > 0:
                    log("SKIP", f"{{sample}}/{{f.name}} already filtered")
                    continue
                if MIN_LEN > 0:
                    kept, dropped = filter_fasta(f, out, MIN_LEN)
                    total_kept += kept
                    total_dropped += dropped
                    if kept == 0:
                        log("WARN", f"{{sample}}/{{f.name}} kept 0 contigs (>= {{MIN_LEN}} bp)")
                else:
                    shutil.copy2(f, out)
            log("OK", f"{{sample}}: {{len(members)}} assemblies organized")

        log("OK", f"contig filter: kept={{total_kept}} dropped={{total_dropped}} (threshold {{MIN_LEN}} bp)")
        PYEOF

        log_warn_errors_from "$STAGE_LOG"
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s08 compress --------------------------------------------------------

def stage_compress(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s08_compress"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    # Only emit --max_contigs if non-default; preserves compatibility with
    # older autocycler versions that don't accept the flag.
    extra = ""
    if int(cfg["max_contigs"]) != 25:
        extra = f" --max_contigs {int(cfg['max_contigs'])}"
    for s in samples:
        sample_dir = f"{layout.autocycler}/{s.sample}"
        # autocycler compress writes input_assemblies.fasta into the autocycler dir
        marker = f"{sample_dir}/input_assemblies.fasta"
        body += textwrap.dedent(
            f"""
            if [ -s {marker} ]; then
                log_event SKIP "{s.sample} (input_assemblies.fasta exists)"
            else
                set +e
                autocycler compress \\
                    --threads {cfg['sge_pe_parallel']} \\
                    --assemblies_dir {layout.assemblies_fasta}/{s.sample} \\
                    --autocycler_dir {sample_dir}{extra}
                rc=$?
                set -e
                if [ $rc -ne 0 ]; then
                    log_event WARN "{s.sample}: compress rc=$rc"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    log_event OK "{s.sample} compressed"
                fi
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s09 cluster ---------------------------------------------------------

def stage_cluster(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s09_cluster"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += f"mkdir -p {layout.cluster_dir}\n\n"
    for s in samples:
        ac_dir = f"{layout.autocycler}/{s.sample}"
        # autocycler cluster writes clustering/ under the autocycler dir
        marker = f"{ac_dir}/clustering"
        body += textwrap.dedent(
            f"""
            if [ -d {marker} ]; then
                log_event SKIP "{s.sample} (clustering/ exists)"
            elif [ ! -s {ac_dir}/input_assemblies.fasta ]; then
                log_event WARN "{s.sample}: compress output missing, skipping"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                set +e
                autocycler cluster --autocycler_dir {ac_dir}
                rc=$?
                set -e
                if [ $rc -ne 0 ]; then
                    log_event WARN "{s.sample}: cluster rc=$rc"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    log_event OK "{s.sample} clustered"
                fi
            fi
            mkdir -p {layout.cluster_dir}/{s.sample}
            ln -sfn {ac_dir}/clustering {layout.cluster_dir}/{s.sample}/clustering

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s10 trim + resolve --------------------------------------------------

def stage_trim_resolve(cfg: dict, layout: Layout, hold: str) -> str:
    name = "s10_trim_resolve"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += textwrap.dedent(
        f"""
        find {layout.cluster_dir} -type d -path "*/clustering/qc_pass" | while read qc_dir; do
            for c in "$qc_dir"/cluster_*; do
                [ -d "$c" ] || continue
                if [ -f "$c/2_trimmed.gfa" ] || [ -f "$c/3_trimmed.gfa" ] || [ -f "$c/trimmed.gfa" ]; then
                    log_event SKIP "trim $c (already trimmed)"
                else
                    set +e
                    autocycler trim -c "$c"
                    rc=$?
                    set -e
                    if [ $rc -ne 0 ]; then
                        log_event WARN "trim $c rc=$rc"
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    else
                        log_event OK "trim $c"
                    fi
                fi
            done
        done

        find {layout.cluster_dir} -type d -path "*/clustering/qc_pass" | while read qc_dir; do
            for c in "$qc_dir"/cluster_*; do
                [ -d "$c" ] || continue
                if [ -f "$c/5_final.gfa" ]; then
                    log_event SKIP "resolve $c (5_final.gfa exists)"
                else
                    set +e
                    autocycler resolve -c "$c"
                    rc=$?
                    set -e
                    if [ $rc -ne 0 ]; then
                        log_event WARN "resolve $c rc=$rc"
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    else
                        log_event OK "resolve $c"
                    fi
                fi
            done
        done

        log_warn_errors_from "$STAGE_LOG"
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s11 combine ---------------------------------------------------------

def stage_combine(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s11_combine"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_autocycler"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    for s in samples:
        ac_dir = f"{layout.autocycler}/{s.sample}"
        marker = f"{ac_dir}/consensus_assembly.fasta"
        body += textwrap.dedent(
            f"""
            if [ -s {marker} ]; then
                log_event SKIP "{s.sample} (consensus_assembly.fasta exists)"
            else
                gfas=$(ls {ac_dir}/clustering/qc_pass/cluster_*/5_final.gfa 2>/dev/null || true)
                if [ -z "$gfas" ]; then
                    log_event WARN "{s.sample}: no 5_final.gfa, skipping combine"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    set +e
                    autocycler combine -a {ac_dir} -i $gfas
                    rc=$?
                    set -e
                    if [ $rc -ne 0 ]; then
                        log_event WARN "{s.sample}: combine rc=$rc"
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    else
                        log_event OK "{s.sample} combined"
                    fi
                fi
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s12 collect ---------------------------------------------------------

def stage_collect(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s12_collect"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        1, cfg["email"], hold_jid=hold,
    )
    body += log_setup(name, layout) + "\n\n"
    body += f"mkdir -p {layout.consensus}\n\n"
    for s in samples:
        dst = f"{layout.consensus}/{s.sample}.fasta"
        body += textwrap.dedent(
            f"""
            if [ -s {dst} ]; then
                log_event SKIP "{s.sample} (already collected)"
            else
                src={layout.autocycler}/{s.sample}/consensus_assembly.fasta
                if [ ! -f "$src" ]; then
                    src={layout.cluster_dir}/{s.sample}/consensus_assembly.fasta
                fi
                if [ -f "$src" ]; then
                    cp "$src" {dst}
                    log_event OK "{s.sample} collected from $src"
                else
                    log_event WARN "{s.sample}: no consensus_assembly.fasta found"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                fi
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s13 annotate --------------------------------------------------------

def stage_annotate(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s13_annotate"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_prokka"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += f"mkdir -p {layout.prokka}\n\n"
    for s in samples:
        proteins = f"--proteins {s.reference_proteins}" if s.reference_proteins else ""
        out_dir = f"{layout.prokka}/{s.sample}"
        marker = f"{out_dir}/{s.sample}.gff"
        body += textwrap.dedent(
            f"""
            if [ -s {marker} ]; then
                log_event SKIP "{s.sample} (prokka .gff exists)"
            elif [ ! -s {layout.consensus}/{s.sample}.fasta ]; then
                log_event WARN "{s.sample}: consensus fasta missing"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                # prokka refuses to overwrite an existing --outdir
                rm -rf {out_dir}
                set +e
                prokka --dbdir {cfg['prokka_db_dir']} \\
                    --outdir {out_dir} \\
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
                rc=$?
                set -e
                if [ $rc -ne 0 ]; then
                    log_event WARN "{s.sample}: prokka rc=$rc"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    log_event OK "{s.sample} annotated"
                fi
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    body += '\n# Final pipeline summary\n'
    body += textwrap.dedent(
        f"""
        log_event SUMMARY "pipeline finished — log at $PIPELINE_LOG"
        echo "" >> $PIPELINE_LOG
        echo "==================== pipeline summary ====================" >> $PIPELINE_LOG
        grep -E '\\[(START|DONE|FAIL)\\]' $PIPELINE_LOG | tail -200 >> $PIPELINE_LOG
        echo "==========================================================" >> $PIPELINE_LOG
        """
    )
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def write_submit_all(layout: Layout, ordered_jobs: List[List[str]]) -> None:
    lines = [
        "#!/bin/bash",
        "# Auto-generated by generate_pipeline.py",
        "# Submits every stage with -hold_jid wired to the previous stage.",
        "# Safe to re-run: each stage checks for existing outputs and skips work",
        "# that's already done. To force a re-run of a sample, delete its outputs.",
        "set -euo pipefail",
        f"cd {layout.submit_dir}",
        "",
        f"# Append a 'pipeline started' marker to the unified log.",
        f"echo \"$(date -u +%Y-%m-%dT%H:%M:%SZ) [START] submit_all $(whoami)@$(hostname)\" >> {layout.root}/pipeline.log",
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
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--pod5-dir", required=True, type=Path)
    p.add_argument("--samplesheet", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--config", required=True, type=Path)
    args = p.parse_args()

    if not args.pod5_dir.is_dir():
        raise SystemExit(f"--pod5-dir does not exist: {args.pod5_dir}")
    if not args.samplesheet.is_file():
        raise SystemExit(f"--samplesheet not found: {args.samplesheet}")
    if not args.config.is_file():
        raise SystemExit(f"--config not found: {args.config}")

    cfg = load_config(args.config)
    samples = parse_samplesheet(args.samplesheet)
    layout = Layout.from_root(args.output_dir.resolve(), args.pod5_dir.resolve())

    for d in [
        layout.root, layout.submit_dir, layout.stdout_dir,
        layout.basecalled, layout.fastq_raw, layout.fastq,
        layout.fastq_filtered, layout.autocycler, layout.assemblies_flat,
        layout.assemblies_fasta, layout.cluster_dir, layout.consensus,
        layout.prokka,
    ]:
        d.mkdir(parents=True, exist_ok=True)

    # Touch the pipeline log so tail -f works immediately
    (layout.root / "pipeline.log").touch()

    n01 = stage_basecall(cfg, layout)
    n02 = stage_demux(cfg, layout, hold=n01)
    n03 = stage_rename(cfg, layout, samples, hold=n02)
    n04 = stage_filter(cfg, layout, samples, hold=n03)
    n05 = stage_subsample(cfg, layout, samples, hold=n04)
    n06_list = stage_assemble(cfg, layout, samples, hold=n05)
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

    print(f"Generated pipeline at: {layout.root}")
    print(f"  - submit_scripts/   ({sum(len(s) for s in ordered)} scripts)")
    print(f"  - submit_all.sh     (driver)")
    print(f"  - pipeline.log      (live, one line per event)")
    print(f"  - {len(samples)} sample(s): {[s.sample for s in samples]}")
    print()
    print(f"Run with:")
    print(f"  bash {layout.root}/submit_all.sh")
    print()
    print(f"Watch progress with:")
    print(f"  tail -f {layout.root}/pipeline.log")


if __name__ == "__main__":
    main()
