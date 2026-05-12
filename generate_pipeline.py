#!/usr/bin/env python3
"""
generate_pipeline.py — emit an end-to-end Nanopore assembly + annotation
pipeline as a stack of SGE submit scripts wired together with -hold_jid.

Inputs
------
  --pod5-dir      Folder of *.pod5 files (raw signal data)
  --samplesheet   TSV with required columns: barcode, sample, genome_length,
                  genus, species, strain. Optional columns (any may be
                  absent): kingdom, reference_proteins, gram (pos/neg/unknown),
                  plasmid (name).
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
  s13_annotate       Bakta per sample using samplesheet metadata
  s14_classify       GTDB-Tk classify_wf across all Bakta .fna outputs

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
]
OPTIONAL_COLS = [
    "kingdom",             # not used by bakta; kept for samplesheet back-compat
    "reference_proteins",  # bakta --proteins <path>
    "gram",                # bakta --gram + / - / ?  (pos / neg / unknown)
    "plasmid",             # bakta --plasmid <name>
]
GRAM_MAP = {"pos": "+", "neg": "-", "unknown": "?"}

# Bakta requires plasmid names to start with a lowercase 'p' and be INSDC-
# compliant (alphanumeric + dot/dash/underscore, max length ~24). This catches
# the most common mistake (using a sample name as the plasmid name) at parse
# time instead of mid-run.
import re as _re
PLASMID_RE = _re.compile(r"^p[a-zA-Z0-9._-]{0,23}$")


@dataclass
class Sample:
    barcode: str
    sample: str
    genome_length: int
    genus: str
    species: str
    strain: str
    kingdom: str = ""
    reference_proteins: str = ""  # empty -> no --proteins flag
    gram: str = ""                # empty -> no --gram flag
    plasmid: str = ""             # empty -> no --plasmid flag


def parse_samplesheet(path: Path) -> List[Sample]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        missing = [c for c in REQUIRED_COLS if c not in cols]
        if missing:
            raise SystemExit(
                f"Samplesheet {path} missing required columns: {missing}\n"
                f"Required: {REQUIRED_COLS}\n"
                f"Optional: {OPTIONAL_COLS}"
            )
        samples: List[Sample] = []
        for row_num, row in enumerate(reader, start=2):
            present_cols = REQUIRED_COLS + [c for c in OPTIONAL_COLS if c in cols]
            if not any((row.get(c) or "").strip() for c in present_cols):
                continue
            try:
                gl = int((row["genome_length"] or "").replace(",", "").strip())
            except ValueError:
                raise SystemExit(
                    f"Row {row_num}: genome_length is not an integer "
                    f"({row['genome_length']!r})"
                )
            gram = (row.get("gram") or "").strip().lower()
            if gram and gram not in GRAM_MAP:
                raise SystemExit(
                    f"Row {row_num}: gram must be one of {list(GRAM_MAP)} "
                    f"(got {gram!r})"
                )
            plasmid = (row.get("plasmid") or "").strip()
            if plasmid and not PLASMID_RE.match(plasmid):
                raise SystemExit(
                    f"Row {row_num}: plasmid name {plasmid!r} is not Bakta-compatible.\n"
                    f"Bakta requires the name to start with a lowercase 'p' and "
                    f"contain only alphanumerics, dots, dashes, or underscores "
                    f"(max 24 chars). Example: pTp114_TS_ori, pUC19, pV3_KS2.\n"
                    f"Leave the column blank for samples that don't have a plasmid."
                )
            samples.append(
                Sample(
                    barcode=row["barcode"].strip(),
                    sample=row["sample"].strip(),
                    genome_length=gl,
                    genus=row["genus"].strip(),
                    species=row["species"].strip(),
                    strain=row["strain"].strip(),
                    kingdom=(row.get("kingdom") or "").strip(),
                    reference_proteins=(row.get("reference_proteins") or "").strip(),
                    gram=gram,
                    plasmid=plasmid,
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
            # Only match definitive error markers. We drop Traceback because
            # Python-based assemblers (nextDenovo especially) emit tracebacks
            # during normal retry behaviour — real failures still surface as
            # non-zero return codes, which the per-command wrapper logs as WARN.
            # Negative filter excludes:
            #   - canu algorithm names containing "error" (errorRate, error
            #     detection, error adjustment, error threshold)
            #   - nextDenovo retry chatter (nextDenovo, cns_align, Re-write
            #     workdir, empty job, log.critical, .tmpfrunt, sgs_fofn,
            #     hifi_fofn, Delete task:)
            #   - plassembler INFO/WARNING noise about canu/unicycler/no plasmids
            #   - our own log events ([SCAN], [WARN], etc.) to prevent
            #     feedback loops when stage stderr is captured into stage log
            grep -nE '(ERROR:|Error:|FATAL:|Fatal:|FAIL:|panic:|Segmentation fault|core dumped|^Killed$)' "$f" 2>/dev/null \\
                | grep -viE '(error[_ ]?rate|error[_ ]detection|error[_ ]adjustment|error[_ ]threshold|fraction[_ ]error|error[_ ]correction|nextDenovo|cns_align|Re-write workdir|empty job|log\\.critical|\\.tmpfrunt|sgs_fofn|hifi_fofn|Delete task:|plassembler|Unicycler has failed|no plasmids|uncorrected reads|Canu failed to correct|\\[(SCAN|WARN|OK|SKIP|START|DONE|FAIL|INFO|SUMMARY)\\])' \\
                | head -20 \\
                | sed -E 's/\\x1B\\[[0-9;]*[mGK]//g' \\
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
    annotate: Path
    mapping: Path
    meta_flye: Path
    bins: Path
    checkm2: Path
    gtdbtk: Path
    reports: Path
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
            annotate=root / "bakta",
            mapping=root / "mapping",
            meta_flye=root / "meta_flye",
            bins=root / "bins",
            checkm2=root / "checkm2",
            gtdbtk=root / "gtdbtk",
            reports=root / "reports",
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
                if [ -s {out_prefix}.fasta ] || [ -f {out_prefix}.done ]; then
                    log_event SKIP "{label} (already complete)"
                else
                    log_event OK "{label} starting"
                    set +e
                    {cmd}
                    rc=$?
                    set -e
                    if [ $rc -ne 0 ]; then
                        log_event WARN "{label} assembler exited with code $rc"
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    elif [ ! -s {out_prefix}.fasta ]; then
                        # rc=0 + no output usually means the assembler ran fine
                        # but found nothing to assemble — typical for plassembler
                        # on samples with no plasmids. Not a failure. Touch a
                        # .done sentinel so future restarts don't re-run it.
                        touch {out_prefix}.done
                        log_event INFO "{label} completed with no output (e.g. plassembler with no plasmids)"
                    else
                        # Also drop a .done marker for samples that did produce
                        # output, so the SKIP check on resume is uniform.
                        touch {out_prefix}.done
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
    max_contigs = int(cfg["max_contigs"])
    body += textwrap.dedent(
        f"""
        log_event OK "organizing assemblies (min_contig_length={min_len}, max_contigs_per_assembly={max_contigs})"
        python3 - <<'PYEOF'
        import sys, shutil
        from collections import defaultdict
        from pathlib import Path

        src = Path("{layout.assemblies_flat}")
        dst = Path("{layout.assemblies_fasta}")
        dst.mkdir(parents=True, exist_ok=True)
        MIN_LEN = {min_len}
        MAX_CONTIGS = {max_contigs}
        PIPELINE_LOG = Path("{layout.root}/pipeline.log")

        def log(level, msg):
            import datetime
            try:
                ts = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
            except AttributeError:
                ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            line = f"{{ts}} [{{level}}] s07_organize {{msg}}\\n"
            with PIPELINE_LOG.open("a") as fh:
                fh.write(line)
            sys.stderr.write(line)

        ASSEMBLERS = {{"canu","flye","metamdbg","miniasm","necat",
                       "nextdenovo","plassembler","raven"}}

        def read_fasta(path):
            \"\"\"Yield (header, sequence) tuples from a fasta file.\"\"\"
            header = None
            seq = []
            with path.open() as fh:
                for line in fh:
                    line = line.rstrip("\\n")
                    if line.startswith(">"):
                        if header is not None:
                            yield header, "".join(seq)
                        header = line
                        seq = []
                    else:
                        seq.append(line)
            if header is not None:
                yield header, "".join(seq)

        def count_contigs(path):
            n = 0
            with path.open() as fh:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
            return n

        def filter_and_cap(in_path, out_path, min_len, max_contigs):
            \"\"\"Drop contigs shorter than min_len, then keep only the top
            max_contigs longest remaining contigs. Writes to out_path.\"\"\"
            kept_short_filter = []
            dropped_short = 0
            for header, seq in read_fasta(in_path):
                if len(seq) >= min_len:
                    kept_short_filter.append((header, seq))
                else:
                    dropped_short += 1
            # Sort descending by length, take top max_contigs
            kept_short_filter.sort(key=lambda x: len(x[1]), reverse=True)
            capped_out = max_contigs > 0 and len(kept_short_filter) > max_contigs
            final = kept_short_filter[:max_contigs] if max_contigs > 0 else kept_short_filter
            dropped_cap = len(kept_short_filter) - len(final)
            with out_path.open("w") as fout:
                for header, seq in final:
                    fout.write(header + "\\n")
                    for i in range(0, len(seq), 80):
                        fout.write(seq[i:i+80] + "\\n")
            return len(final), dropped_short, dropped_cap, capped_out

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
        total_dropped_short = 0
        total_dropped_cap = 0
        capped_count = 0
        for sample, members in groups.items():
            sub = dst / sample
            sub.mkdir(exist_ok=True)
            for f in members:
                out = sub / f.name
                if out.exists() and out.stat().st_size > 0:
                    # Skip only if existing output already satisfies the cap.
                    # If it has more contigs than the cap, re-filter to apply
                    # the new threshold.
                    existing = count_contigs(out)
                    if MAX_CONTIGS <= 0 or existing <= MAX_CONTIGS:
                        log("SKIP", f"{{sample}}/{{f.name}} already filtered ({{existing}} contigs)")
                        continue
                    else:
                        log("INFO", f"{{sample}}/{{f.name}} has {{existing}} contigs > cap {{MAX_CONTIGS}}, refiltering")
                kept, ds, dc, capped = filter_and_cap(f, out, MIN_LEN, MAX_CONTIGS)
                total_kept += kept
                total_dropped_short += ds
                total_dropped_cap += dc
                if capped:
                    capped_count += 1
                if kept == 0:
                    log("WARN", f"{{sample}}/{{f.name}} kept 0 contigs after filter (min_len={{MIN_LEN}})")
            log("OK", f"{{sample}}: {{len(members)}} assemblies organized")

        log("OK", f"contig filter: kept={{total_kept}} dropped_short={{total_dropped_short}} dropped_over_cap={{total_dropped_cap}} (capped {{capped_count}} assemblies at {{MAX_CONTIGS}})")
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
        # autocycler compress writes input_assemblies.gfa into the autocycler dir
        marker = f"{sample_dir}/input_assemblies.gfa"
        body += textwrap.dedent(
            f"""
            if [ -s {marker} ]; then
                log_event SKIP "{s.sample} (input_assemblies.gfa exists)"
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
            elif [ ! -s {ac_dir}/input_assemblies.gfa ]; then
                log_event WARN "{s.sample}: compress output (input_assemblies.gfa) missing, skipping"
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
        # Search the real autocycler dir, not the cluster/ symlink. find
        # without -L doesn't traverse symlinks, so searching through
        # cluster/<sample>/clustering (a symlink) would never see qc_pass.
        find {layout.autocycler} -type d -path "*/clustering/qc_pass" -not -path "*/cluster/*" | while read qc_dir; do
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

        find {layout.autocycler} -type d -path "*/clustering/qc_pass" -not -path "*/cluster/*" | while read qc_dir; do
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
    name = "s16_annotate"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_bakta"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += f"mkdir -p {layout.annotate}\n\n"
    bakta_min_len = int(cfg.get("bakta_min_contig_length", 200))
    bakta_threads = int(cfg.get("bakta_threads", 20))
    for s in samples:
        # Build optional flags conditionally so they're only present when
        # the samplesheet supplied a value. Each flag is on its own line with
        # the same continuation indentation as the other bakta flags so the
        # generated script reads cleanly.
        opt_flags: List[str] = []
        if s.plasmid:
            opt_flags.append(f"--plasmid {s.plasmid}")
        if s.gram:
            opt_flags.append(f"--gram {GRAM_MAP[s.gram]}")
        if s.reference_proteins:
            opt_flags.append(f"--proteins {s.reference_proteins}")
        # Indented to match the other --flag lines after textwrap.dedent
        # strips the common leading whitespace.
        opt_lines = "".join(f"        {flag} \\\n" for flag in opt_flags)

        out_dir = f"{layout.annotate}/{s.sample}"
        marker = f"{out_dir}/{s.sample}.gff3"
        body += textwrap.dedent(
            f"""\
            if [ -s {marker} ]; then
                log_event SKIP "{s.sample} (bakta .gff3 exists)"
            elif [ ! -s {layout.consensus}/{s.sample}.fasta ]; then
                log_event WARN "{s.sample}: consensus fasta missing"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                set +e
                bakta \\
                    --db {cfg['bakta_db']} \\
                    --min-contig-length {bakta_min_len} \\
                    --prefix {s.sample} \\
                    --output {out_dir} \\
                    --genus {s.genus} \\
                    --species {s.species} \\
                    --strain {s.strain} \\
            """
        ) + opt_lines + textwrap.dedent(
            f"""\
                    --threads {bakta_threads} \\
                    --force \\
                    {layout.consensus}/{s.sample}.fasta
                rc=$?
                set -e
                if [ $rc -ne 0 ]; then
                    log_event WARN "{s.sample}: bakta rc=$rc"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    log_event OK "{s.sample} annotated"
                fi
            fi

            """
        )
    # Also run Bakta on every MetaBAT2 bin. Bins are putative MAGs of
    # unknown taxonomy, so we pass only --strain <bin_id> for traceability
    # and let Bakta annotate without --genus/--species/--proteins/--gram.
    # Output goes to bakta/<bin_id>/ in parallel with the consensus outputs.
    body += textwrap.dedent(
        f"""

        # ---- Bakta on MetaBAT2 bins ----
        for bin_fa in {layout.bins}/*/*.fa; do
            [ -e "$bin_fa" ] || continue
            bin_id=$(basename "$bin_fa" .fa)
            out_dir={layout.annotate}/${{bin_id}}
            marker=${{out_dir}}/${{bin_id}}.gff3
            if [ -s "$marker" ]; then
                log_event SKIP "${{bin_id}} (bakta .gff3 exists)"
            else
                set +e
                bakta \\
                    --db {cfg['bakta_db']} \\
                    --min-contig-length {bakta_min_len} \\
                    --prefix "$bin_id" \\
                    --output "$out_dir" \\
                    --strain "$bin_id" \\
                    --threads {bakta_threads} \\
                    --force \\
                    "$bin_fa"
                rc=$?
                set -e
                if [ $rc -ne 0 ]; then
                    log_event WARN "$bin_id: bakta rc=$rc"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    log_event OK "$bin_id annotated (bin)"
                fi
            fi
        done

        """
    )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s13 map reads (minimap2 + samtools) ---------------------------------

def stage_map_reads(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s13_map_reads"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_mapping"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    threads = int(cfg.get("minimap2_threads", 20))
    body += f"mkdir -p {layout.mapping}\n\n"
    for s in samples:
        sdir = f"{layout.mapping}/{s.sample}"
        marker = f"{sdir}/flagstat.txt"
        body += textwrap.dedent(
            f"""
            mkdir -p {sdir}
            if [ -s {marker} ] && [ -s {sdir}/unmapped.fastq ]; then
                log_event SKIP "{s.sample} (already mapped)"
            elif [ ! -s {layout.consensus}/{s.sample}.fasta ]; then
                log_event WARN "{s.sample}: consensus fasta missing, skipping"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            elif [ ! -s {layout.fastq_filtered}/{s.sample}.fastq ]; then
                log_event WARN "{s.sample}: filtered reads missing, skipping"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                set +e
                minimap2 -ax map-ont -t {threads} \\
                    {layout.consensus}/{s.sample}.fasta \\
                    {layout.fastq_filtered}/{s.sample}.fastq 2>/dev/null \\
                    | samtools sort -@ {threads} -o {sdir}/mapped.bam -
                rc=$?
                if [ $rc -eq 0 ]; then
                    samtools index {sdir}/mapped.bam
                    samtools flagstat {sdir}/mapped.bam > {marker}
                    # Extract unmapped reads as fastq (-f 4 selects unmapped)
                    samtools fastq -f 4 {sdir}/mapped.bam > {sdir}/unmapped.fastq 2>/dev/null
                    mapped=$(grep -E '[0-9]+ \\+ [0-9]+ mapped \\(' {marker} | head -1 | awk '{{print $1}}')
                    total=$(head -1 {marker} | awk '{{print $1}}')
                    unmapped=$(grep -cE '^@|^>' {sdir}/unmapped.fastq 2>/dev/null || echo 0)
                    unmapped_reads=$(awk 'NR%4==1' {sdir}/unmapped.fastq | wc -l)
                    log_event OK "{s.sample} mapped=$mapped/$total unmapped_reads=$unmapped_reads"
                else
                    log_event WARN "{s.sample}: minimap2/samtools rc=$rc"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                fi
                set -e
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s15 meta-assemble unmapped reads (Flye --meta) ---------------------

def stage_meta_assemble(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s14_meta_assemble"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_assembly"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_flye_meta"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    threads = int(cfg.get("flye_meta_threads", 20))
    min_reads = int(cfg.get("flye_meta_min_unmapped_reads", 100))
    body += f"mkdir -p {layout.meta_flye}\n\n"
    for s in samples:
        unmapped = f"{layout.mapping}/{s.sample}/unmapped.fastq"
        out_dir = f"{layout.meta_flye}/{s.sample}"
        marker = f"{out_dir}/assembly.fasta"
        skip_sentinel = f"{out_dir}/.skipped_low_reads"
        body += textwrap.dedent(
            f"""
            if [ -s {marker} ]; then
                log_event SKIP "{s.sample} (meta assembly exists)"
            elif [ -f {skip_sentinel} ]; then
                log_event SKIP "{s.sample} (previously skipped — too few unmapped reads)"
            elif [ ! -s {unmapped} ]; then
                log_event INFO "{s.sample}: no unmapped reads, skipping meta-assembly"
                mkdir -p {out_dir} && touch {skip_sentinel}
            else
                n_reads=$(awk 'NR%4==1' {unmapped} | wc -l)
                if [ "$n_reads" -lt {min_reads} ]; then
                    log_event INFO "{s.sample}: only $n_reads unmapped reads (< {min_reads}), skipping meta-assembly"
                    mkdir -p {out_dir} && touch {skip_sentinel}
                else
                    log_event OK "{s.sample}: starting flye --meta on $n_reads unmapped reads"
                    rm -rf {out_dir}
                    set +e
                    flye --meta --nano-hq {unmapped} \\
                        --out-dir {out_dir} \\
                        --threads {threads}
                    rc=$?
                    set -e
                    if [ $rc -ne 0 ] || [ ! -s {marker} ]; then
                        log_event WARN "{s.sample}: flye --meta rc=$rc (assembly may be empty)"
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                        mkdir -p {out_dir} && touch {skip_sentinel}
                    else
                        contigs=$(grep -c '^>' {marker})
                        log_event OK "{s.sample}: flye --meta produced $contigs contigs"
                    fi
                fi
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s16 MetaBAT2 binning -----------------------------------------------

def stage_bin(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s15_bin"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_metabat"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    threads = int(cfg.get("minimap2_threads", 20))
    min_contig = int(cfg.get("metabat_min_contig", 1500))
    body += f"mkdir -p {layout.bins}\n\n"
    for s in samples:
        meta_asm = f"{layout.meta_flye}/{s.sample}/assembly.fasta"
        bin_dir = f"{layout.bins}/{s.sample}"
        depth_bam = f"{bin_dir}/depth.bam"
        depth_txt = f"{bin_dir}/depth.txt"
        marker_done = f"{bin_dir}/.binning_done"
        body += textwrap.dedent(
            f"""
            mkdir -p {bin_dir}
            if [ -f {marker_done} ]; then
                log_event SKIP "{s.sample} (binning previously completed)"
            elif [ ! -s {meta_asm} ]; then
                log_event INFO "{s.sample}: no meta-assembly, skipping binning"
                touch {marker_done}
            else
                log_event OK "{s.sample}: mapping unmapped reads to meta-assembly for depth"
                set +e
                minimap2 -ax map-ont -t {threads} \\
                    {meta_asm} \\
                    {layout.mapping}/{s.sample}/unmapped.fastq 2>/dev/null \\
                    | samtools sort -@ {threads} -o {depth_bam} -
                samtools index {depth_bam}
                jgi_summarize_bam_contig_depths --outputDepth {depth_txt} {depth_bam}
                metabat2 \\
                    -i {meta_asm} \\
                    -a {depth_txt} \\
                    -o {bin_dir}/{s.sample}_bin \\
                    -m {min_contig} \\
                    -t {threads}
                rc=$?
                set -e
                # metabat2 may produce 0 bins (returns 0 in that case)
                n_bins=$(ls {bin_dir}/{s.sample}_bin.*.fa 2>/dev/null | wc -l)
                if [ $rc -ne 0 ]; then
                    log_event WARN "{s.sample}: metabat2 rc=$rc"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    log_event OK "{s.sample}: metabat2 produced $n_bins bins"
                fi
                touch {marker_done}
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s17 CheckM2 quality assessment -------------------------------------

def stage_checkm2(cfg: dict, layout: Layout, hold: str) -> str:
    name = "s17_checkm2"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_gtdbtk"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_checkm2"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    threads = int(cfg.get("checkm2_threads", 20))
    db = cfg.get("checkm2_db", "") or ""
    out_dir = f"{layout.checkm2}/output"
    input_dir = f"{layout.checkm2}/input"
    marker = f"{out_dir}/quality_report.tsv"
    if db:
        body += f'export CHECKM2DB="{db}"\n'
    body += textwrap.dedent(
        f"""

        if [ -s {marker} ]; then
            log_event SKIP "checkm2 quality_report.tsv already exists"
        else
            # Stage all assemblies (consensuses + bins) into a flat dir with .fasta extension
            mkdir -p {input_dir}
            rm -f {input_dir}/*.fasta
            # Consensus assemblies from autocycler
            for f in {layout.consensus}/*.fasta; do
                [ -e "$f" ] || continue
                ln -sf "$f" {input_dir}/$(basename "$f")
            done
            # MetaBAT2 bins (already include sample name in filename)
            for f in {layout.bins}/*/*.fa; do
                [ -e "$f" ] || continue
                base=$(basename "$f" .fa)
                ln -sf "$f" {input_dir}/${{base}}.fasta
            done
            n_inputs=$(ls {input_dir}/*.fasta 2>/dev/null | wc -l)
            log_event OK "CheckM2 input: $n_inputs assemblies"
            if [ "$n_inputs" -eq 0 ]; then
                log_event WARN "no assemblies found for CheckM2"
                exit 0
            fi
            mkdir -p {out_dir}
            set +e
            checkm2 predict \\
                --threads {threads} \\
                --input {input_dir} \\
                --output-directory {out_dir} \\
                --force \\
                --extension fasta
            rc=$?
            set -e
            if [ $rc -ne 0 ]; then
                log_event WARN "checkm2 predict rc=$rc"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                n_done=$(tail -n +2 {marker} 2>/dev/null | wc -l)
                log_event OK "CheckM2 finished: $n_done assemblies assessed"
            fi
        fi

        log_warn_errors_from "$STAGE_LOG"
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s18 GTDB-Tk classification (consensuses + bins) --------------------

def stage_classify(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s18_classify"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"],
        cfg.get("sge_h_vmem_gtdbtk", "100G"),
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_gtdbtk"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"

    ext = cfg.get("gtdbtk_extension", "fna")
    cpus = int(cfg.get("gtdbtk_cpus", 20))
    pplacer_cpus = int(cfg.get("gtdbtk_pplacer_cpus", 12))
    data_path = cfg.get("gtdbtk_data_path", "") or ""

    out_dir = f"{layout.gtdbtk}/output"
    batchfile = f"{layout.gtdbtk}/batchfile.tsv"

    body += f"mkdir -p {layout.gtdbtk}\n"
    if data_path:
        body += f'export GTDBTK_DATA_PATH="{data_path}"\n'
    body += textwrap.dedent(
        f"""

        # Build the GTDB-Tk batchfile from BOTH Bakta annotated consensuses AND
        # MetaBAT2 bins. Format: <fasta_path>\\t<genome_id>
        # Bin IDs are <sample>_bin.N so origin sample is identifiable downstream.
        : > {batchfile}
        """
    )
    for s in samples:
        fna = f"{layout.annotate}/{s.sample}/{s.sample}.{ext}"
        body += textwrap.dedent(
            f"""
            if [ -s {fna} ]; then
                printf '%s\\t%s\\n' '{fna}' '{s.sample}' >> {batchfile}
            else
                log_event WARN "{s.sample}: missing bakta output ({fna}), excluded from GTDB-Tk"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            fi
            """
        )
    body += textwrap.dedent(
        f"""

        # Add every MetaBAT2 bin. Prefer the Bakta-annotated .{ext} (cleaner
        # headers), fall back to the raw .fa if Bakta hasn't run on the bin.
        # Genome ID is <sample>_bin.N so origin is traceable in GTDB-Tk rows.
        for f in {layout.bins}/*/*.fa; do
            [ -e "$f" ] || continue
            genome_id=$(basename "$f" .fa)
            bakta_fna={layout.annotate}/${{genome_id}}/${{genome_id}}.{ext}
            if [ -s "$bakta_fna" ]; then
                printf '%s\\t%s\\n' "$bakta_fna" "$genome_id" >> {batchfile}
            else
                printf '%s\\t%s\\n' "$f" "$genome_id" >> {batchfile}
            fi
        done

        n_entries=$(wc -l < {batchfile})
        log_event OK "GTDB-Tk batchfile has $n_entries entries (consensuses + bins)"
        if [ "$n_entries" -eq 0 ]; then
            log_event WARN "no Bakta outputs found, skipping GTDB-Tk entirely"
            exit 0
        fi

        # Resume check: GTDB-Tk writes gtdbtk.bac120.summary.tsv and/or
        # gtdbtk.ar53.summary.tsv at the very end. If either exists, skip.
        if compgen -G "{out_dir}/gtdbtk.*.summary.tsv" > /dev/null; then
            log_event SKIP "GTDB-Tk summary already present in {out_dir}"
        else
            mkdir -p {out_dir}
            set +e
            gtdbtk classify_wf \\
                --batchfile {batchfile} \\
                --out_dir {out_dir} \\
                --extension {ext} \\
                --cpus {cpus} \\
                --pplacer_cpus {pplacer_cpus}
            rc=$?
            set -e
            if [ $rc -ne 0 ]; then
                log_event WARN "gtdbtk classify_wf rc=$rc"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                summary_files=$(ls {out_dir}/gtdbtk.*.summary.tsv 2>/dev/null)
                log_event OK "GTDB-Tk classification complete (summaries: $summary_files)"
            fi
        fi

        log_warn_errors_from "$STAGE_LOG"
        """
    )
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s19 aggregate metrics ----------------------------------------------

def stage_aggregate(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s19_aggregate"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        1, cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_report"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += f"mkdir -p {layout.reports}\n\n"
    sample_names = [s.sample for s in samples]
    body += textwrap.dedent(
        f"""
        log_event OK "collecting metrics across stages"
        python3 - <<'PYEOF'
        import csv, json, re, sys
        from pathlib import Path

        ROOT = Path("{layout.root}")
        OUT  = Path("{layout.reports}") / "metrics.tsv"
        SAMPLES = {sample_names!r}

        def fastq_count(p):
            if not p.exists() or p.stat().st_size == 0:
                return 0
            n = 0
            with p.open() as fh:
                for i, _ in enumerate(fh):
                    if i % 4 == 0:
                        n += 1
            return n

        def fasta_stats(p):
            if not p.exists() or p.stat().st_size == 0:
                return 0, 0, 0
            contigs, total, longest = 0, 0, 0
            cur = 0
            with p.open() as fh:
                for line in fh:
                    if line.startswith(">"):
                        if cur > 0:
                            total += cur
                            longest = max(longest, cur)
                            contigs += 1
                        cur = 0
                    else:
                        cur += len(line.strip())
                if cur > 0:
                    total += cur
                    longest = max(longest, cur)
                    contigs += 1
            return contigs, total, longest

        def fasta_per_contig(p):
            \"\"\"Return list of (header_name, length) tuples, sorted descending by length.
            For autocycler consensuses, the largest contig is the chromosome and
            the rest are plasmids/extrachromosomal elements.\"\"\"
            if not p.exists() or p.stat().st_size == 0:
                return []
            out = []
            name = None
            cur = 0
            with p.open() as fh:
                for line in fh:
                    if line.startswith(">"):
                        if name is not None:
                            out.append((name, cur))
                        name = line[1:].split()[0] if line[1:].strip() else "unnamed"
                        cur = 0
                    else:
                        cur += len(line.strip())
                if name is not None:
                    out.append((name, cur))
            out.sort(key=lambda x: x[1], reverse=True)
            return out

        def parse_flagstat(p):
            if not p.exists():
                return None
            txt = p.read_text()
            total = mapped = 0
            m = re.match(r"^(\\d+)", txt)
            if m: total = int(m.group(1))
            m = re.search(r"^(\\d+)\\s*\\+\\s*\\d+\\s+mapped\\s*\\(", txt, re.M)
            if m: mapped = int(m.group(1))
            return total, mapped

        def load_bakta_json(p):
            try:
                with p.open() as fh:
                    return json.load(fh)
            except Exception:
                return {{}}

        def bakta_feature_counts(p):
            \"\"\"Try Bakta JSON first; fall back to .txt summary if needed.\"\"\"
            sample_dir = p.parent
            stem = p.stem
            j = sample_dir / f"{{stem}}.json"
            counts = {{"CDS": 0, "tRNA": 0, "rRNA": 0, "tmRNA": 0,
                      "ncRNA": 0, "CRISPR": 0, "sORF": 0, "oriC": 0, "oriT": 0}}
            if j.exists():
                d = load_bakta_json(j)
                feats = d.get("features", [])
                for f in feats:
                    t = f.get("type", "")
                    if t in counts:
                        counts[t] += 1
                    elif t == "ncRNA-region":
                        counts["ncRNA"] += 1
                return counts
            # Fall back to .txt summary if no JSON
            t = sample_dir / f"{{stem}}.txt"
            if t.exists():
                for line in t.read_text().splitlines():
                    for k in counts:
                        if line.startswith(k + ":") or line.startswith(k + " "):
                            try:
                                counts[k] = int(line.split()[-1])
                            except Exception:
                                pass
            return counts

        def load_checkm2():
            p = ROOT / "checkm2" / "output" / "quality_report.tsv"
            if not p.exists():
                return {{}}
            d = {{}}
            with p.open() as fh:
                rdr = csv.DictReader(fh, delimiter="\\t")
                for r in rdr:
                    d[r["Name"]] = {{
                        "completeness": float(r.get("Completeness", 0) or 0),
                        "contamination": float(r.get("Contamination", 0) or 0),
                    }}
            return d

        def load_gtdbtk():
            d = {{}}
            for f in (ROOT / "gtdbtk" / "output").glob("gtdbtk.*.summary.tsv"):
                with f.open() as fh:
                    rdr = csv.DictReader(fh, delimiter="\\t")
                    for r in rdr:
                        gid = r.get("user_genome", "")
                        d[gid] = {{
                            "classification": r.get("classification", ""),
                            "domain":  r.get("classification", "").split(";")[0] if ";" in r.get("classification","") else "",
                        }}
            return d

        def parse_taxonomy(classification):
            \"\"\"Parse a GTDB-Tk classification string into ranked components.\"\"\"
            levels = {{"d__": "domain", "p__": "phylum", "c__": "class",
                       "o__": "order", "f__": "family", "g__": "genus",
                       "s__": "species"}}
            out = {{v: "" for v in levels.values()}}
            for p in (classification or "").split(";"):
                p = p.strip()
                for prefix, name in levels.items():
                    if p.startswith(prefix):
                        out[name] = p[len(prefix):]
                        break
            # GTDB species fields look like "Escherichia coli" — drop genus prefix
            if out["species"]:
                out["species_short"] = out["species"]
            else:
                out["species_short"] = ""
            return out

        def load_expected():
            p = Path("{layout.reports}") / "expected_taxonomy.tsv"
            d = {{}}
            if p.exists():
                with p.open() as fh:
                    rdr = csv.DictReader(fh, delimiter="\\t")
                    for r in rdr:
                        d[r["sample"]] = {{
                            "expected_genus":   r.get("expected_genus", ""),
                            "expected_species": r.get("expected_species", ""),
                            "expected_strain":  r.get("expected_strain", ""),
                        }}
            return d

        checkm2  = load_checkm2()
        gtdbtk   = load_gtdbtk()
        expected = load_expected()

        rows = []
        for s in SAMPLES:
            raw_reads  = fastq_count(ROOT / "fastq" / f"{{s}}.fastq")
            filt_reads = fastq_count(ROOT / "fastq_filtered" / f"{{s}}.fastq")
            flag = parse_flagstat(ROOT / "mapping" / s / "flagstat.txt")
            total, mapped = (flag if flag else (0, 0))
            unmapped_reads = fastq_count(ROOT / "mapping" / s / "unmapped.fastq")
            cons = ROOT / "consensus" / f"{{s}}.fasta"
            c_n, c_len, c_long = fasta_stats(cons)
            # Split consensus into chromosome (largest contig) + plasmids (rest)
            contigs_sorted = fasta_per_contig(cons)
            if contigs_sorted:
                chrom_name, chrom_len = contigs_sorted[0]
                plasmid_contigs = contigs_sorted[1:]
            else:
                chrom_name, chrom_len = "", 0
                plasmid_contigs = []
            plasmid_count = len(plasmid_contigs)
            plasmid_total_bp = sum(L for _, L in plasmid_contigs)
            plasmid_names = ";".join(n for n, _ in plasmid_contigs)
            plasmid_lengths = ";".join(str(L) for _, L in plasmid_contigs)
            bakta_fna = ROOT / "bakta" / s / f"{{s}}.fna"
            feats = bakta_feature_counts(bakta_fna) if bakta_fna.exists() else {{}}
            meta_asm = ROOT / "meta_flye" / s / "assembly.fasta"
            m_n, m_len, _ = fasta_stats(meta_asm)
            bin_files = sorted((ROOT / "bins" / s).glob(f"{{s}}_bin.*.fa")) if (ROOT / "bins" / s).exists() else []
            n_bins = len(bin_files)
            cm = checkm2.get(s, {{}})
            gt = gtdbtk.get(s, {{}})
            exp = expected.get(s, {{}})
            tax = parse_taxonomy(gt.get("classification", ""))
            exp_g = exp.get("expected_genus", "")
            exp_sp = exp.get("expected_species", "")
            # Genus match: GTDB's genus == samplesheet's genus
            genus_match = (tax["genus"] != "" and exp_g != ""
                           and tax["genus"].lower() == exp_g.lower())
            # Species match: GTDB's species is "Genus species" — compare full string
            exp_full = f"{{exp_g}} {{exp_sp}}".strip()
            species_match = (tax["species"] != "" and exp_full != ""
                             and tax["species"].lower() == exp_full.lower())

            base = {{
                "sample": s,
                "assembly_type": "isolate_consensus",
                "genome_id": s,
                "raw_reads": raw_reads,
                "filtered_reads": filt_reads,
                "total_reads_mapped_input": total,
                "mapped_reads": mapped,
                "unmapped_reads": unmapped_reads,
                "prop_mapped": (mapped / total if total else 0),
                "prop_unmapped": (1 - mapped / total) if total else 0,
                "n_contigs": c_n,
                "genome_length_bp": c_len,
                "longest_contig_bp": c_long,
                "chromosome_length_bp": chrom_len,
                "chromosome_contig": chrom_name,
                "plasmid_count": plasmid_count,
                "plasmid_total_bp": plasmid_total_bp,
                "plasmid_names": plasmid_names,
                "plasmid_lengths": plasmid_lengths,
                "meta_assembly_contigs": m_n,
                "meta_assembly_length_bp": m_len,
                "n_bins": n_bins,
                "completeness": cm.get("completeness", ""),
                "contamination": cm.get("contamination", ""),
                "classification": gt.get("classification", ""),
                "domain": gt.get("domain", ""),
                "phylum": tax["phylum"],
                "class": tax["class"],
                "order": tax["order"],
                "family": tax["family"],
                "genus": tax["genus"],
                "species": tax["species"],
                "expected_genus": exp_g,
                "expected_species": exp_sp,
                "genus_match": genus_match,
                "species_match": species_match,
                "cds": feats.get("CDS", 0),
                "tRNA": feats.get("tRNA", 0),
                "rRNA": feats.get("rRNA", 0),
                "tmRNA": feats.get("tmRNA", 0),
                "ncRNA": feats.get("ncRNA", 0),
                "CRISPR": feats.get("CRISPR", 0),
            }}
            rows.append(base)

            # One row per bin
            for bf in bin_files:
                bid = bf.stem  # e.g. <sample>_bin.1
                bn, blen, blong = fasta_stats(bf)
                cmb = checkm2.get(bid, {{}})
                gtb = gtdbtk.get(bid, {{}})
                btax = parse_taxonomy(gtb.get("classification", ""))
                rows.append({{
                    "sample": s,
                    "assembly_type": "metagenomic_bin",
                    "genome_id": bid,
                    "n_contigs": bn,
                    "genome_length_bp": blen,
                    "longest_contig_bp": blong,
                    "completeness": cmb.get("completeness", ""),
                    "contamination": cmb.get("contamination", ""),
                    "classification": gtb.get("classification", ""),
                    "domain": gtb.get("domain", ""),
                    "phylum": btax["phylum"],
                    "class": btax["class"],
                    "order": btax["order"],
                    "family": btax["family"],
                    "genus": btax["genus"],
                    "species": btax["species"],
                    "expected_genus": exp_g,
                    "expected_species": exp_sp,
                }})

        fields = ["sample","assembly_type","genome_id","raw_reads",
                  "filtered_reads","total_reads_mapped_input","mapped_reads",
                  "unmapped_reads","prop_mapped","prop_unmapped","n_contigs",
                  "genome_length_bp","longest_contig_bp",
                  "chromosome_length_bp","chromosome_contig",
                  "plasmid_count","plasmid_total_bp","plasmid_names",
                  "plasmid_lengths","meta_assembly_contigs",
                  "meta_assembly_length_bp","n_bins","completeness",
                  "contamination","classification","domain","phylum","class",
                  "order","family","genus","species","expected_genus",
                  "expected_species","genus_match","species_match",
                  "cds","tRNA","rRNA","tmRNA","ncRNA","CRISPR"]
        with OUT.open("w") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\\t",
                               extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"Wrote {{len(rows)}} rows to {{OUT}}")
        PYEOF

        n_rows=$(tail -n +2 {layout.reports}/metrics.tsv | wc -l)
        log_event OK "aggregated $n_rows rows into metrics.tsv"
        log_warn_errors_from "$STAGE_LOG"
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# Standalone Python script for the report stage. Written to submit_scripts/
# so that the bash wrapper can just invoke it without an embedded heredoc
# (and so the report code stays out of the outer-f-string brace minefield).
REPORT_PY = r'''#!/usr/bin/env python3
"""Generate publication-ready plots + HTML report from metrics.tsv."""
import argparse, base64, sys
from html import escape
from pathlib import Path

try:
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError as e:
    sys.stderr.write(f"missing dependency: {e}\n")
    sys.exit(1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metrics", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--plots-dir", required=True, type=Path)
    args = p.parse_args()

    if not args.metrics.exists():
        sys.stderr.write("metrics.tsv not found — run s19_aggregate first\n")
        sys.exit(1)

    df = pd.read_csv(args.metrics, sep="\t")
    if df.empty:
        sys.stderr.write("metrics.tsv is empty\n")
        sys.exit(1)

    args.plots_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "figure.figsize": (10, 5.5),
        "figure.dpi": 130,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.0,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "legend.frameon": False,
        "font.family": "DejaVu Sans",
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": "--",
    })

    COL_ISOLATE = "#2c7fb8"
    COL_BIN     = "#fdae61"
    COL_MAPPED  = "#41ab5d"
    COL_UNMAP   = "#d73027"

    def save_plot(name):
        path = args.plots_dir / f"{name}.png"
        plt.savefig(path)
        plt.close()
        return path

    isolates = df[df["assembly_type"] == "isolate_consensus"].copy()
    bins     = df[df["assembly_type"] == "metagenomic_bin"].copy()
    plots = []

    # 1. Read flow per sample
    if not isolates.empty:
        x = np.arange(len(isolates))
        w = 0.2
        fig, ax = plt.subplots()
        ax.bar(x - 1.5*w, isolates["raw_reads"],      w, label="Raw reads",      color="#bdbdbd")
        ax.bar(x - 0.5*w, isolates["filtered_reads"], w, label="After chopper",  color="#737373")
        ax.bar(x + 0.5*w, isolates["mapped_reads"],   w, label="Mapped",         color=COL_MAPPED)
        ax.bar(x + 1.5*w, isolates["unmapped_reads"], w, label="Unmapped",       color=COL_UNMAP)
        ax.set_xticks(x)
        ax.set_xticklabels(isolates["sample"], rotation=45, ha="right")
        ax.set_ylabel("Read count")
        ax.set_title("Read flow per sample")
        ax.legend(loc="upper right", ncol=2)
        plots.append(("Read flow per sample",
                      "Raw demultiplexed reads, post-chopper QC reads, reads mapping "
                      "to the Autocycler consensus, and reads remaining unmapped "
                      "(input to metagenomic Flye).",
                      save_plot("read_flow")))

    # 2. Mapping proportion
    if not isolates.empty:
        fig, ax = plt.subplots()
        x = np.arange(len(isolates))
        pm = pd.to_numeric(isolates["prop_mapped"], errors="coerce").fillna(0).values
        pu = 1 - pm
        ax.bar(x, pm, color=COL_MAPPED, label="Mapped")
        ax.bar(x, pu, bottom=pm, color=COL_UNMAP, label="Unmapped")
        ax.set_xticks(x)
        ax.set_xticklabels(isolates["sample"], rotation=45, ha="right")
        ax.set_ylabel("Proportion of reads")
        ax.set_title("Mapping proportion per sample")
        ax.set_ylim(0, 1)
        ax.legend(loc="upper right")
        plots.append(("Mapping proportion",
                      "Fraction of QC-passing reads aligning to the per-sample "
                      "consensus. Unmapped fraction is fed into metagenomic Flye "
                      "+ MetaBAT2 to recover co-assembled bins.",
                      save_plot("prop_mapped")))

    # 3. Assembly sizes
    if not isolates.empty:
        fig, ax = plt.subplots()
        x = np.arange(len(isolates))
        cons_len = pd.to_numeric(isolates["genome_length_bp"], errors="coerce").fillna(0) / 1e6
        sample_bin_len = (
            bins.groupby("sample")["genome_length_bp"]
                .apply(lambda s: pd.to_numeric(s, errors="coerce").sum())
            if not bins.empty else pd.Series(dtype=float)
        ) / 1e6
        bin_vals = [sample_bin_len.get(s, 0) for s in isolates["sample"]]
        ax.bar(x, cons_len, color=COL_ISOLATE, label="Isolate consensus")
        ax.bar(x, bin_vals, bottom=cons_len, color=COL_BIN, label="MetaBAT2 bins (total)")
        ax.set_xticks(x)
        ax.set_xticklabels(isolates["sample"], rotation=45, ha="right")
        ax.set_ylabel("Genome size (Mb)")
        ax.set_title("Assembly size per sample")
        ax.legend(loc="upper right")
        plots.append(("Assembly size per sample",
                      "Length of Autocycler consensus plus total length of MetaBAT2 "
                      "bins recovered from unmapped reads.",
                      save_plot("assembly_size")))

    # 4. CheckM2 completeness vs contamination
    qc_cols = ["genome_id", "assembly_type", "completeness", "contamination"]
    qc = pd.concat([
        isolates[qc_cols] if all(c in isolates.columns for c in qc_cols) else pd.DataFrame(columns=qc_cols),
        bins[qc_cols]     if all(c in bins.columns     for c in qc_cols) else pd.DataFrame(columns=qc_cols),
    ])
    qc["completeness"]  = pd.to_numeric(qc["completeness"],  errors="coerce")
    qc["contamination"] = pd.to_numeric(qc["contamination"], errors="coerce")
    qc = qc.dropna(subset=["completeness","contamination"])
    if not qc.empty:
        fig, ax = plt.subplots()
        for typ, col, mk in [("isolate_consensus", COL_ISOLATE, "o"),
                             ("metagenomic_bin",   COL_BIN,     "s")]:
            d = qc[qc["assembly_type"] == typ]
            if not d.empty:
                ax.scatter(d["completeness"], d["contamination"],
                           c=col, marker=mk, s=70, edgecolor="white",
                           linewidth=0.8, label=typ.replace("_", " "),
                           alpha=0.85)
        ax.axhline(5,  color="grey", linestyle=":", linewidth=0.8)
        ax.axvline(90, color="grey", linestyle=":", linewidth=0.8)
        ax.text(91, max(qc["contamination"].max()*0.92, 4.5),
                "MIMAG high-quality\n(>=90% / <=5%)", fontsize=8, color="grey")
        ax.set_xlabel("Completeness (%)")
        ax.set_ylabel("Contamination (%)")
        ax.set_title("CheckM2 genome quality")
        ax.set_xlim(0, 105)
        ax.set_ylim(-0.5, max(15, qc["contamination"].max() + 2))
        ax.legend(loc="upper left")
        plots.append(("CheckM2 quality",
                      "Completeness vs contamination for every assembly. "
                      "Dashed lines mark MIMAG high-quality thresholds.",
                      save_plot("checkm2")))

    # 5. Annotation feature counts
    if not isolates.empty:
        feat_cols = ["cds","tRNA","rRNA","tmRNA","ncRNA","CRISPR"]
        present = [c for c in feat_cols if c in isolates.columns]
        if present:
            fig, ax = plt.subplots()
            x = np.arange(len(isolates))
            width = 0.8 / len(present)
            palette = ["#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a","#b15928"]
            any_drawn = False
            for i, c in enumerate(present):
                vals = pd.to_numeric(isolates[c], errors="coerce").fillna(0)
                if vals.sum() == 0:
                    continue
                ax.bar(x + i*width - 0.4 + width/2, vals, width,
                       label=c, color=palette[i % len(palette)])
                any_drawn = True
            if any_drawn:
                ax.set_xticks(x)
                ax.set_xticklabels(isolates["sample"], rotation=45, ha="right")
                ax.set_ylabel("Feature count")
                ax.set_title("Bakta annotation features per sample")
                ax.legend(loc="upper right", ncol=3)
                plots.append(("Bakta annotation features",
                              "Coding sequences and non-coding RNA features called "
                              "by Bakta in each isolate consensus.",
                              save_plot("annotation_features")))
            else:
                plt.close()

    # 6. CDS density (CDS vs genome size)
    if not isolates.empty:
        try:
            cds = pd.to_numeric(isolates["cds"], errors="coerce")
            gl  = pd.to_numeric(isolates["genome_length_bp"], errors="coerce") / 1e6
            ok  = cds.notna() & gl.notna() & (cds > 0) & (gl > 0)
            if ok.sum() > 1:
                fig, ax = plt.subplots()
                gl_ok = gl[ok].astype(float).values
                cds_ok = cds[ok].astype(float).values
                sample_ok = isolates["sample"][ok].values
                ax.scatter(gl_ok, cds_ok, c=COL_ISOLATE, s=70,
                           edgecolor="white", linewidth=0.8)
                for sample, x, y in zip(sample_ok, gl_ok, cds_ok):
                    ax.annotate(str(sample), (float(x), float(y)),
                                fontsize=7, xytext=(3,3), textcoords="offset points",
                                color="#555")
                ax.set_xlabel("Genome size (Mb)")
                ax.set_ylabel("CDS count")
                ax.set_title("Coding density (CDS vs genome size)")
                plots.append(("Coding density",
                              "Bakta CDS count vs genome size. Bacteria typically "
                              "fall close to a ~900 CDS/Mb line.",
                              save_plot("cds_vs_size")))
        except Exception as e:
            sys.stderr.write(f"plot cds_vs_size skipped: {e}\n")
            plt.close("all")

    # 7. Taxonomy at phylum level
    tax = pd.concat([
        isolates[["genome_id","assembly_type","classification"]] if "classification" in isolates.columns else pd.DataFrame(),
        bins[["genome_id","assembly_type","classification"]]     if "classification" in bins.columns     else pd.DataFrame(),
    ])
    tax = tax[tax["classification"].fillna("") != ""]
    if not tax.empty:
        def phylum(c):
            for p in (c or "").split(";"):
                if p.startswith("p__"):
                    return p.replace("p__", "") or "Unclassified"
            return "Unclassified"
        tax = tax.assign(phylum=tax["classification"].apply(phylum))
        counts = tax.groupby(["phylum","assembly_type"]).size().unstack(fill_value=0)
        counts = counts.loc[counts.sum(axis=1).sort_values(ascending=False).index]
        fig, ax = plt.subplots(figsize=(10, max(4, 0.4*len(counts))))
        counts.plot(kind="barh", stacked=True, ax=ax,
                    color=[COL_ISOLATE, COL_BIN][:counts.shape[1]])
        ax.set_xlabel("Number of assemblies")
        ax.set_title("GTDB-Tk classification (phylum level)")
        ax.invert_yaxis()
        ax.legend(loc="lower right")
        plots.append(("Taxonomy (phylum)",
                      "GTDB-Tk phylum-level assignments across all assemblies.",
                      save_plot("taxonomy_phylum")))

    # 8. Bins per sample
    if not bins.empty:
        counts = bins.groupby("sample").size()
        fig, ax = plt.subplots()
        ax.bar(np.arange(len(counts)), counts.values, color=COL_BIN)
        ax.set_xticks(np.arange(len(counts)))
        ax.set_xticklabels(counts.index, rotation=45, ha="right")
        ax.set_ylabel("Number of MetaBAT2 bins")
        ax.set_title("Metagenomic bins per sample")
        plots.append(("Bins per sample",
                      "Number of putative MAGs recovered by MetaBAT2 from the "
                      "Flye --meta assembly of unmapped reads.",
                      save_plot("bins_per_sample")))

    # ============================================================
    # Comparative plots: expected vs observed taxonomy
    # ============================================================

    all_asm = pd.concat([isolates, bins], ignore_index=True) if not isolates.empty or not bins.empty else pd.DataFrame()
    for c in ["genus", "species", "phylum", "expected_genus", "expected_species"]:
        if c not in all_asm.columns:
            all_asm[c] = ""
    all_asm[["genus","species","phylum","expected_genus","expected_species"]] = \
        all_asm[["genus","species","phylum","expected_genus","expected_species"]].fillna("").astype(str)

    # 9. Expected vs observed — sample-level genus/species match heat-table
    if not isolates.empty:
        try:
            iso_view = isolates[["sample","expected_genus","expected_species","genus","species","genus_match","species_match"]].copy()
            iso_view = iso_view.astype(object).fillna("")
            fig_h = max(3.5, 0.45 * len(iso_view) + 1)
            fig, ax = plt.subplots(figsize=(11, fig_h))
            ax.set_axis_off()
            col_labels = ["Sample", "Expected genus", "Expected species",
                          "Observed genus (GTDB)", "Observed species (GTDB)", "Genus match"]
            cell_text = []
            cell_colors = []
            for _, r in iso_view.iterrows():
                gm = str(r.get("genus_match","")).lower() in ("true","1")
                sm = str(r.get("species_match","")).lower() in ("true","1")
                exp_sp = (r["expected_genus"] + " " + r["expected_species"]).strip()
                row = [r["sample"], r["expected_genus"] or "—",
                       exp_sp or "—",
                       r["genus"] or "—",
                       r["species"] or "—",
                       ("species" if sm else ("genus" if gm else "mismatch"))]
                cell_text.append(row)
                if sm:
                    row_color = "#c7e9b4"   # green: species match
                elif gm:
                    row_color = "#ffffb2"   # yellow: genus match only
                else:
                    row_color = "#fcae91"   # red: mismatch
                cell_colors.append(["#ffffff", "#ffffff", "#ffffff",
                                    "#ffffff", "#ffffff", row_color])
            tbl = ax.table(cellText=cell_text, colLabels=col_labels,
                           cellColours=cell_colors, loc="center", cellLoc="left",
                           colWidths=[0.16,0.17,0.20,0.17,0.20,0.10])
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(8.5)
            tbl.scale(1, 1.35)
            for c in range(len(col_labels)):
                tbl[(0, c)].set_facecolor("#1a3a5c")
                tbl[(0, c)].set_text_props(color="white", weight="bold")
            ax.set_title("Expected vs observed taxonomy (isolate consensuses)",
                         fontsize=12, fontweight="bold", pad=10)
            plots.append(("Expected vs observed taxonomy",
                          "Per-sample comparison of the samplesheet's expected "
                          "genus/species against the GTDB-Tk classification of the "
                          "Autocycler consensus. Green = species match, yellow = "
                          "genus match only, red = mismatch.",
                          save_plot("expected_vs_observed")))
        except Exception as e:
            sys.stderr.write(f"plot expected_vs_observed skipped: {e}\n")
            plt.close("all")

    # 10. Per-sample species composition — stacked bar of assembled Mb by species
    if not all_asm.empty:
        try:
            tmp = all_asm.copy()
            tmp["genome_length_bp"] = pd.to_numeric(tmp["genome_length_bp"], errors="coerce").fillna(0)
            def label(row):
                sp = row.get("species","")
                if sp:
                    return sp
                gn = row.get("genus","")
                if gn:
                    return f"{gn} sp."
                ph = row.get("phylum","")
                return f"Unclassified {ph}" if ph else "Unclassified"
            tmp["taxon_label"] = tmp.apply(label, axis=1)
            comp = tmp.pivot_table(index="sample", columns="taxon_label",
                                   values="genome_length_bp", aggfunc="sum",
                                   fill_value=0) / 1e6
            if comp.shape[1] > 0:
                # Order columns: place expected match first per sample, then others by total size
                col_order = comp.sum(axis=0).sort_values(ascending=False).index.tolist()
                comp = comp[col_order]
                n_taxa = comp.shape[1]
                cmap = plt.cm.tab20
                colors = [cmap(i % 20) for i in range(n_taxa)]
                fig, ax = plt.subplots(figsize=(11, max(5, 0.4*len(comp)+2)))
                comp.plot(kind="barh", stacked=True, ax=ax, color=colors, width=0.78)
                ax.invert_yaxis()
                ax.set_xlabel("Assembled length (Mb)")
                ax.set_title("Per-sample assembly composition by GTDB classification")
                ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0),
                          fontsize=7.5, frameon=False, ncol=1,
                          title="Taxon (most → least abundant)",
                          title_fontsize=8.5)
                plots.append(("Per-sample assembly composition",
                              "Total assembled length per sample, broken down by "
                              "GTDB-Tk classification across all assemblies "
                              "(consensus + bins). Shows the proportional makeup "
                              "of unique taxa recovered from each sample.",
                              save_plot("composition_by_species")))
        except Exception as e:
            sys.stderr.write(f"plot composition_by_species skipped: {e}\n")
            plt.close("all")

    # 11. Sample × phylum heatmap (count of assemblies)
    if not all_asm.empty:
        try:
            tmp = all_asm.copy()
            tmp["phylum_label"] = tmp["phylum"].replace("", "Unclassified")
            mat = tmp.pivot_table(index="sample", columns="phylum_label",
                                  values="genome_id", aggfunc="count",
                                  fill_value=0)
            if mat.size > 0:
                mat = mat[mat.sum(axis=0).sort_values(ascending=False).index]
                fig, ax = plt.subplots(figsize=(max(7, 0.55*mat.shape[1]+3),
                                                max(4, 0.4*mat.shape[0]+1)))
                im = ax.imshow(mat.values, cmap="YlGnBu", aspect="auto")
                ax.set_xticks(np.arange(mat.shape[1]))
                ax.set_xticklabels(mat.columns, rotation=40, ha="right", fontsize=9)
                ax.set_yticks(np.arange(mat.shape[0]))
                ax.set_yticklabels(mat.index, fontsize=9)
                # Annotate cells
                for i in range(mat.shape[0]):
                    for j in range(mat.shape[1]):
                        v = mat.values[i, j]
                        if v > 0:
                            ax.text(j, i, str(int(v)), ha="center", va="center",
                                    fontsize=8.5,
                                    color="white" if v > mat.values.max()*0.55 else "#1a3a5c")
                cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
                cb.set_label("Assemblies", fontsize=9)
                ax.set_title("Phylum × sample assembly count")
                ax.set_xlabel("")
                ax.set_ylabel("")
                plots.append(("Phylum × sample heatmap",
                              "Number of assemblies (consensus + bins) per "
                              "sample-phylum combination. Cells show how the "
                              "assembled biomass distributed across higher taxa "
                              "for each input sample.",
                              save_plot("phylum_heatmap")))
        except Exception as e:
            sys.stderr.write(f"plot phylum_heatmap skipped: {e}\n")
            plt.close("all")

    # 12. Species richness per sample
    if not all_asm.empty:
        try:
            tmp = all_asm.copy()
            tmp = tmp[tmp["species"] != ""]
            if not tmp.empty:
                rich = tmp.groupby("sample")["species"].nunique().sort_values(ascending=False)
                fig, ax = plt.subplots(figsize=(10, max(3.5, 0.35*len(rich)+2)))
                bars = ax.barh(np.arange(len(rich)), rich.values,
                               color=COL_ISOLATE, edgecolor="white", linewidth=0.6)
                ax.set_yticks(np.arange(len(rich)))
                ax.set_yticklabels(rich.index)
                ax.invert_yaxis()
                ax.set_xlabel("Unique GTDB species")
                ax.set_title("Species richness per sample")
                for i, v in enumerate(rich.values):
                    ax.text(v + 0.05, i, str(int(v)), va="center", fontsize=9, color="#444")
                plots.append(("Species richness",
                              "Number of distinct GTDB-Tk species classifications "
                              "recovered per sample, pooling the isolate consensus "
                              "and all metagenomic bins.",
                              save_plot("species_richness")))
        except Exception as e:
            sys.stderr.write(f"plot species_richness skipped: {e}\n")
            plt.close("all")

    # 13. Sample × genus bubble plot
    if not all_asm.empty:
        try:
            tmp = all_asm.copy()
            tmp["genus_label"] = tmp["genus"].replace("", "Unclassified")
            tmp["genome_length_bp"] = pd.to_numeric(tmp["genome_length_bp"], errors="coerce").fillna(0)
            tmp = tmp[tmp["genus_label"] != "Unclassified"]
            if not tmp.empty:
                bub = tmp.groupby(["sample","genus_label"]).agg(
                    n=("genome_id","count"),
                    bp=("genome_length_bp","sum")).reset_index()
                genera = sorted(bub["genus_label"].unique(),
                                key=lambda g: bub[bub["genus_label"] == g]["bp"].sum(),
                                reverse=True)
                samples_ord = sorted(bub["sample"].unique())
                if len(genera) > 0 and len(samples_ord) > 0:
                    fig, ax = plt.subplots(figsize=(max(7, 0.6*len(genera)+3),
                                                    max(4, 0.45*len(samples_ord)+1.5)))
                    max_bp = bub["bp"].max() if not bub.empty else 1
                    for _, r in bub.iterrows():
                        x = genera.index(r["genus_label"])
                        y = samples_ord.index(r["sample"])
                        size = 60 + 700 * (r["bp"] / max_bp)
                        ax.scatter(x, y, s=size, color=COL_ISOLATE,
                                   edgecolor="white", linewidth=1.0, alpha=0.85)
                        ax.text(x, y, str(int(r["n"])), ha="center", va="center",
                                fontsize=8, color="white", fontweight="bold")
                    ax.set_xticks(np.arange(len(genera)))
                    ax.set_xticklabels(genera, rotation=40, ha="right", fontsize=9, style="italic")
                    ax.set_yticks(np.arange(len(samples_ord)))
                    ax.set_yticklabels(samples_ord, fontsize=9)
                    ax.set_xlim(-0.5, len(genera) - 0.5)
                    ax.set_ylim(-0.5, len(samples_ord) - 0.5)
                    ax.invert_yaxis()
                    ax.set_title("Sample × genus distribution (bubble size = assembled bp)")
                    plots.append(("Sample × genus bubble plot",
                                  "Each bubble represents one or more assemblies of a "
                                  "given genus in a given sample. Bubble size scales with "
                                  "total assembled length; numbers inside are assembly "
                                  "counts.",
                                  save_plot("sample_genus_bubble")))
        except Exception as e:
            sys.stderr.write(f"plot sample_genus_bubble skipped: {e}\n")
            plt.close("all")

    # 14. Replicon breakdown — chromosome vs plasmids vs bins per sample
    if not isolates.empty:
        try:
            x = np.arange(len(isolates))
            chrom = pd.to_numeric(isolates.get("chromosome_length_bp", 0), errors="coerce").fillna(0) / 1e6
            plas  = pd.to_numeric(isolates.get("plasmid_total_bp", 0),    errors="coerce").fillna(0) / 1e6
            # Sum bin lengths per sample
            if not bins.empty:
                sample_bin_len = bins.groupby("sample")["genome_length_bp"].apply(
                    lambda s: pd.to_numeric(s, errors="coerce").sum()) / 1e6
                bin_vals = np.array([sample_bin_len.get(s, 0) for s in isolates["sample"]])
            else:
                bin_vals = np.zeros(len(isolates))
            fig, ax = plt.subplots(figsize=(11, max(5, 0.5*len(isolates)+1)))
            ax.barh(x, chrom, color="#2c7fb8", label="Chromosome")
            ax.barh(x, plas,  left=chrom, color="#7b68ee", label="Plasmid(s)")
            ax.barh(x, bin_vals, left=chrom + plas, color=COL_BIN, label="MetaBAT2 bins")
            ax.set_yticks(x)
            ax.set_yticklabels(isolates["sample"])
            ax.invert_yaxis()
            ax.set_xlabel("Assembled length (Mb)")
            ax.set_title("Per-sample replicon breakdown")
            ax.legend(loc="lower right")
            # Annotate counts on the right of each bar
            for i, (cl, pl, bn) in enumerate(zip(chrom, plas, bin_vals)):
                n_plas = int(isolates["plasmid_count"].iloc[i]) if "plasmid_count" in isolates.columns else 0
                n_bin = int(bins[bins["sample"] == isolates["sample"].iloc[i]].shape[0]) if not bins.empty else 0
                ax.text(cl + pl + bn + 0.05, i,
                        f"chr {cl:.2f} Mb · {n_plas} plasmid(s) · {n_bin} bin(s)",
                        fontsize=8, va="center", color="#444")
            plots.append(("Replicon breakdown per sample",
                          "For each isolate: the chromosome (longest consensus "
                          "contig), the sum of plasmid contigs (other consensus "
                          "contigs), and the total length of MetaBAT2 bins recovered "
                          "from unmapped reads.",
                          save_plot("replicon_breakdown")))
        except Exception as e:
            sys.stderr.write(f"plot replicon_breakdown skipped: {e}\n")
            plt.close("all")

    # ---- HTML ----
    def img_tag(p):
        b64 = base64.b64encode(p.read_bytes()).decode()
        return f'<img class="plot" src="data:image/png;base64,{b64}" alt="{p.stem}">'

    def table_html(d):
        return d.to_html(index=False, classes="metrics", border=0, na_rep="—")

    n_iso  = len(isolates)
    n_bins = len(bins)
    n_samples = isolates["sample"].nunique() if not isolates.empty else 0
    total = n_iso + n_bins
    raw_total = int(pd.to_numeric(isolates["raw_reads"], errors="coerce").fillna(0).sum()) if not isolates.empty else 0

    css = """
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
           background: #f8f9fa; color: #222; margin: 0; padding: 2rem 1rem; }
    .container { max-width: 1100px; margin: 0 auto; }
    h1 { color: #1a3a5c; border-bottom: 3px solid #1a3a5c; padding-bottom: .5rem; }
    h2 { color: #1a3a5c; margin-top: 2.5rem; border-bottom: 1px solid #ccc; padding-bottom: .25rem; }
    .meta { color: #666; font-size: .9rem; margin-bottom: 1.5rem; }
    .plot-card { background: white; border-radius: 8px; padding: 1.5rem;
                 box-shadow: 0 1px 3px rgba(0,0,0,.08); margin-bottom: 1.5rem; }
    .plot-title { font-weight: 600; font-size: 1.1rem; color: #1a3a5c; margin-bottom: .25rem; }
    .plot-caption { color: #555; font-size: .9rem; margin-bottom: 1rem; }
    img.plot { width: 100%; height: auto; display: block; }
    table.metrics { border-collapse: collapse; font-size: .82rem; width: 100%; background: white; }
    table.metrics th, table.metrics td { padding: .35rem .6rem; border-bottom: 1px solid #eee; text-align: left; }
    table.metrics th { background: #1a3a5c; color: white; font-weight: 600; }
    table.metrics tr:hover td { background: #f3f6fa; }
    .summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                    gap: 1rem; margin-bottom: 2rem; }
    .summary-card { background: white; padding: 1rem 1.25rem; border-radius: 8px;
                    box-shadow: 0 1px 3px rgba(0,0,0,.08); }
    .summary-card .value { font-size: 1.8rem; font-weight: 700; color: #1a3a5c; }
    .summary-card .label { color: #666; font-size: .85rem; text-transform: uppercase; letter-spacing: .05em; }
    .sample-card { background: white; border-radius: 8px; padding: 1.5rem 1.75rem;
                   box-shadow: 0 1px 3px rgba(0,0,0,.08); margin-bottom: 1.5rem; }
    .sample-card h3 { margin: 0 0 .25rem 0; color: #1a3a5c; font-size: 1.25rem; }
    .sample-card .sample-meta { color: #666; font-size: .9rem; margin: 0 0 1rem 0; }
    .replicon-section { padding: .65rem 0; border-top: 1px solid #eee; }
    .replicon-section:first-of-type { border-top: none; }
    .replicon-label { font-weight: 600; font-size: .85rem;
                      text-transform: uppercase; letter-spacing: .05em;
                      display: inline-block; padding: .15rem .55rem;
                      border-radius: 4px; margin-bottom: .35rem; }
    .replicon-label.main-asm { background: #d9eaf7; color: #1a3a5c; }
    .replicon-label.plasmids { background: #ece6f7; color: #4b3a87; }
    .replicon-label.bins     { background: #fff1d6; color: #8a5a00; }
    .replicon-list { margin: .25rem 0 .5rem 1.25rem; padding: 0; }
    .replicon-list li { margin: .25rem 0; font-size: .9rem; line-height: 1.4; }
    .replicon-empty { color: #888; font-style: italic; font-size: .9rem; margin: .25rem 0; }
    .bin-class { display: inline-block; margin-top: .15rem; font-size: .82rem; color: #555; }
    .unclass { color: #999; font-style: italic; }
    .badge { display: inline-block; margin-left: .5rem; padding: .12rem .55rem;
             background: #e8f5e9; color: #1b5e20; font-size: .8rem; border-radius: 4px; }
    code { background: #f3f6fa; padding: .05rem .35rem; border-radius: 3px;
           font-size: .85rem; }
    """

    cards = (
        f'<div class="summary-grid">'
        f'<div class="summary-card"><div class="value">{n_samples}</div><div class="label">Samples</div></div>'
        f'<div class="summary-card"><div class="value">{n_iso}</div><div class="label">Isolate consensuses</div></div>'
        f'<div class="summary-card"><div class="value">{n_bins}</div><div class="label">MetaBAT2 bins</div></div>'
        f'<div class="summary-card"><div class="value">{total}</div><div class="label">Total assemblies</div></div>'
        f'<div class="summary-card"><div class="value">{raw_total:,}</div><div class="label">Raw reads (sum)</div></div>'
        f'</div>'
    )

    plot_blocks = []
    for title, caption, png in plots:
        plot_blocks.append(
            f'<div class="plot-card">'
            f'<div class="plot-title">{escape(title)}</div>'
            f'<div class="plot-caption">{escape(caption)}</div>'
            f'{img_tag(png)}'
            f'</div>'
        )

    # ---- Per-sample summary cards (main assembly + plasmids + bins) ----
    def fmt_bp(v):
        try:
            v = float(v)
            if v >= 1e6:
                return f"{v/1e6:.2f} Mb"
            if v >= 1e3:
                return f"{v/1e3:.0f} kb"
            return f"{int(v)} bp"
        except Exception:
            return "—"

    def fmt_int(v, default="—"):
        try:
            return f"{int(float(v)):,}"
        except Exception:
            return default

    def fmt_classification(c):
        if not c:
            return "<span class='unclass'>not classified</span>"
        # Bold genus + species portion if present
        parts = c.split(";")
        last = parts[-1].strip() if parts else ""
        return f"<code>{escape(c)}</code>"

    sample_cards = []
    if not isolates.empty:
        sorted_samples = isolates["sample"].tolist()
        bin_by_sample = bins.groupby("sample") if not bins.empty else None
        for s in sorted_samples:
            row = isolates[isolates["sample"] == s].iloc[0]
            def _safe_str(v):
                if v is None or (isinstance(v, float) and pd.isna(v)):
                    return ""
                return str(v)
            classification = _safe_str(row.get("classification"))
            chrom_bp = row.get("chromosome_length_bp", 0)
            try:
                n_plas = int(float(row.get("plasmid_count", 0) or 0))
            except (TypeError, ValueError):
                n_plas = 0
            plas_names_raw = _safe_str(row.get("plasmid_names"))
            plas_lens_raw  = _safe_str(row.get("plasmid_lengths"))
            plas_names = plas_names_raw.split(";") if plas_names_raw else []
            plas_lens  = plas_lens_raw.split(";")  if plas_lens_raw  else []
            plas_list  = [(n, l) for n, l in zip(plas_names, plas_lens) if n]
            exp_gen    = _safe_str(row.get("expected_genus"))
            exp_sp     = _safe_str(row.get("expected_species"))
            exp_str    = (exp_gen + " " + exp_sp).strip() or "—"
            comp = row.get("completeness", "")
            cont = row.get("contamination", "")
            # Render plasmid block
            if plas_list:
                plas_rows = "".join(
                    f"<li><code>{escape(n)}</code> &mdash; {fmt_bp(l)}</li>"
                    for n, l in plas_list
                )
                plas_html = f"<ul class='replicon-list'>{plas_rows}</ul>"
            elif n_plas > 0:
                plas_html = f"<p class='replicon-empty'>{n_plas} plasmid contig(s) (unnamed)</p>"
            else:
                plas_html = "<p class='replicon-empty'>None detected.</p>"

            # Bins for this sample
            if bin_by_sample is not None and s in bin_by_sample.groups:
                bdf = bin_by_sample.get_group(s)
                bin_rows = []
                for _, b in bdf.iterrows():
                    bgid = b.get("genome_id", "")
                    blen = b.get("genome_length_bp", 0)
                    bclass = b.get("classification", "") or ""
                    bcomp = b.get("completeness", "")
                    bcont = b.get("contamination", "")
                    qc_str = ""
                    if bcomp != "" and pd.notna(bcomp):
                        qc_str = f" &middot; CheckM2 {fmt_int(bcomp, '?')}%/{fmt_int(bcont, '?')}%"
                    bin_rows.append(
                        f"<li><code>{escape(str(bgid))}</code> &mdash; {fmt_bp(blen)}"
                        f"{qc_str}<br><span class='bin-class'>{fmt_classification(bclass)}</span></li>"
                    )
                bin_html = f"<ul class='replicon-list'>{''.join(bin_rows)}</ul>"
            else:
                bin_html = "<p class='replicon-empty'>None recovered.</p>"

            # CheckM2 badge for main assembly
            if comp != "" and pd.notna(comp):
                badge = f"<span class='badge'>CheckM2: {fmt_int(comp,'?')}% complete &middot; {fmt_int(cont,'?')}% contamination</span>"
            else:
                badge = ""

            sample_cards.append(
                f"<div class='sample-card'>"
                f"<h3>{escape(s)}</h3>"
                f"<p class='sample-meta'>Expected: <em>{escape(exp_str)}</em></p>"
                f"<div class='replicon-section'>"
                f"<div class='replicon-label main-asm'>Main assembly (Autocycler consensus)</div>"
                f"<p>Chromosome contig: {fmt_bp(chrom_bp)}{' ' + badge if badge else ''}</p>"
                f"<p>Classification: {fmt_classification(classification)}</p>"
                f"</div>"
                f"<div class='replicon-section'>"
                f"<div class='replicon-label plasmids'>Plasmid contigs ({n_plas})</div>"
                f"{plas_html}"
                f"</div>"
                f"<div class='replicon-section'>"
                f"<div class='replicon-label bins'>MetaBAT2 bins ({int(row.get('n_bins',0) or 0)})</div>"
                f"{bin_html}"
                f"</div>"
                f"</div>"
            )
    sample_cards_html = "".join(sample_cards) if sample_cards else "<p>No samples processed.</p>"

    iso_cols = ["sample","genome_length_bp","n_contigs","mapped_reads",
                "unmapped_reads","completeness","contamination",
                "classification","cds","tRNA","rRNA","n_bins"]
    iso_cols = [c for c in iso_cols if c in isolates.columns]
    iso_disp = isolates[iso_cols].copy() if not isolates.empty else pd.DataFrame()
    if not iso_disp.empty:
        if "genome_length_bp" in iso_disp:
            iso_disp["genome_length_bp"] = pd.to_numeric(iso_disp["genome_length_bp"], errors="coerce").apply(
                lambda v: f"{int(v):,}" if pd.notna(v) else "—")
        for c in ["mapped_reads","unmapped_reads","cds","tRNA","rRNA"]:
            if c in iso_disp:
                iso_disp[c] = pd.to_numeric(iso_disp[c], errors="coerce").apply(
                    lambda v: f"{int(v):,}" if pd.notna(v) else "—")
        for c in ["completeness","contamination"]:
            if c in iso_disp:
                iso_disp[c] = pd.to_numeric(iso_disp[c], errors="coerce").apply(
                    lambda v: f"{v:.1f}" if pd.notna(v) else "—")

    bin_cols = ["genome_id","sample","genome_length_bp","n_contigs",
                "completeness","contamination","classification"]
    bin_cols = [c for c in bin_cols if c in bins.columns]
    bin_disp = bins[bin_cols].copy() if not bins.empty else pd.DataFrame()
    if not bin_disp.empty:
        if "genome_length_bp" in bin_disp:
            bin_disp["genome_length_bp"] = pd.to_numeric(bin_disp["genome_length_bp"], errors="coerce").apply(
                lambda v: f"{int(v):,}" if pd.notna(v) else "—")
        for c in ["completeness","contamination"]:
            if c in bin_disp:
                bin_disp[c] = pd.to_numeric(bin_disp[c], errors="coerce").apply(
                    lambda v: f"{v:.1f}" if pd.notna(v) else "—")

    html = (
        "<!doctype html>"
        "<html><head><meta charset='utf-8'><title>Nanopore pipeline report</title>"
        f"<style>{css}</style></head><body><div class='container'>"
        "<h1>Nanopore assembly + annotation pipeline report</h1>"
        "<p class='meta'>Generated by stage s20_report. "
        "Source data: <code>reports/metrics.tsv</code>.</p>"
        f"{cards}"
        "<h2>Per-sample assemblies</h2>"
        "<p class='meta'>For each sample: the main assembly (Autocycler "
        "consensus) and its GTDB-Tk classification, any plasmid contigs "
        "detected in the consensus, and any MetaBAT2 bins recovered from "
        "unmapped reads (each with its own classification).</p>"
        f"{sample_cards_html}"
        "<h2>Plots</h2>"
        f"{''.join(plot_blocks) if plot_blocks else '<p>No plots generated (no data?).</p>'}"
        "<h2>Tabular summary (main assemblies)</h2>"
        f"{table_html(iso_disp) if not iso_disp.empty else '<p>No data.</p>'}"
        "<h2>Tabular summary (MetaBAT2 bins)</h2>"
        f"{table_html(bin_disp) if not bin_disp.empty else '<p>No bins recovered.</p>'}"
        "</div></body></html>"
    )
    args.out.write_text(html)
    print(f"Wrote {args.out} ({len(plots)} plots)")


if __name__ == "__main__":
    main()
'''


# ---- s20 generate plots + HTML report ----------------------------------

def stage_report(cfg: dict, layout: Layout, hold: str) -> str:
    name = "s20_report"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        1, cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_report"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"

    # Write the standalone report script once (it's self-contained Python,
    # not generated per-run, so we keep the f-string brace minefield out of
    # the bash side entirely).
    py_path = layout.submit_dir / "_report.py"
    py_path.write_text(REPORT_PY)
    py_path.chmod(py_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    body += textwrap.dedent(
        f"""
        log_event OK "generating plots + HTML report"
        set +e
        python3 {py_path} \\
            --metrics {layout.reports}/metrics.tsv \\
            --out {layout.reports}/report.html \\
            --plots-dir {layout.reports}/plots
        rc=$?
        set -e
        if [ $rc -ne 0 ]; then
            log_event WARN "report generator rc=$rc"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
        else
            log_event OK "report written to {layout.reports}/report.html"
        fi

        log_warn_errors_from "$STAGE_LOG"

        log_event SUMMARY "pipeline finished — log at $PIPELINE_LOG, report at {layout.reports}/report.html"
        echo "" >> $PIPELINE_LOG
        echo "==================== pipeline summary ====================" >> $PIPELINE_LOG
        grep -aE '\\[(START|DONE|FAIL)\\]' $PIPELINE_LOG | tail -200 >> $PIPELINE_LOG
        echo "==========================================================" >> $PIPELINE_LOG
        """
    ).strip() + "\n"
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


def write_pipeline_status(layout: Layout) -> None:
    """Standalone helper the user can run anytime to see which stages are
    in-flight, stalled, or done. Diffs START vs DONE events in pipeline.log
    and reports orphans (started but never finished — could mean still
    running OR killed by SGE; SIGKILL bypasses the EXIT trap)."""
    body = textwrap.dedent(
        f"""
        #!/bin/bash
        # pipeline_status.sh — quick check on what each stage is doing.
        #
        # Usage:  bash {layout.root}/pipeline_status.sh
        #
        # Reads pipeline.log and reports:
        #   - stages with [START] but no [DONE]     (running or killed)
        #   - stages with [DONE] but non-zero failures
        #   - per-level counts
        LOG="{layout.root}/pipeline.log"
        if [ ! -f "$LOG" ]; then
            echo "No pipeline.log at $LOG"
            exit 1
        fi

        echo "=== event counts ==="
        grep -aoE '\\[(START|DONE|FAIL|OK|SKIP|WARN|INFO|SCAN|SUMMARY)\\]' "$LOG" \\
            | sort | uniq -c | sort -rn

        echo ""
        echo "=== stages started but not yet DONE ==="
        comm -23 \\
            <(grep -aE '\\[START\\]' "$LOG" | awk '{{print $3}}' | sort -u) \\
            <(grep -aE '\\[DONE\\]'  "$LOG" | awk '{{print $3}}' | sort -u) \\
            | sed 's/^/  /' || true
        if [ $? -ne 0 ] || [ -z "$(comm -23 \\
            <(grep -aE '\\[START\\]' "$LOG" | awk '{{print $3}}' | sort -u) \\
            <(grep -aE '\\[DONE\\]'  "$LOG" | awk '{{print $3}}' | sort -u))" ]; then
            echo "  (none)"
        fi

        echo ""
        echo "=== stages with per-sample failures ==="
        grep -aE '\\[DONE\\].*per-sample failures' "$LOG" || echo "  (none)"

        echo ""
        echo "=== stage-level FAILs ==="
        grep -aE '\\[FAIL\\]' "$LOG" || echo "  (none)"
        """
    ).strip() + "\n"
    out = layout.root / "pipeline_status.sh"
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
        layout.annotate, layout.mapping, layout.meta_flye, layout.bins,
        layout.checkm2, layout.gtdbtk, layout.reports,
    ]:
        d.mkdir(parents=True, exist_ok=True)

    # Touch the pipeline log so tail -f works immediately
    (layout.root / "pipeline.log").touch()

    # Persist the samplesheet's expected taxonomy so the aggregate + report
    # stages can compare expected vs observed (GTDB-Tk) classifications.
    exp_path = layout.reports / "expected_taxonomy.tsv"
    with exp_path.open("w") as fh:
        fh.write("sample\texpected_genus\texpected_species\texpected_strain\texpected_kingdom\n")
        for s in samples:
            fh.write(f"{s.sample}\t{s.genus}\t{s.species}\t{s.strain}\t{s.kingdom}\n")

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
    # Mapping → meta-assembly → binning happen BEFORE Bakta now, so Bakta
    # can annotate both the autocycler consensuses AND the metabat2 bins in
    # a single stage.
    n13 = stage_map_reads(cfg, layout, samples, hold=n12)
    n14 = stage_meta_assemble(cfg, layout, samples, hold=n13)
    n15 = stage_bin(cfg, layout, samples, hold=n14)
    n16 = stage_annotate(cfg, layout, samples, hold=n15)
    n17 = stage_checkm2(cfg, layout, hold=n16)
    n18 = stage_classify(cfg, layout, samples, hold=n17)
    n19 = stage_aggregate(cfg, layout, samples, hold=n18)
    n20 = stage_report(cfg, layout, hold=n19)

    ordered = [
        [n01], [n02], [n03], [n04], [n05],
        n06_list,
        [n07], [n08], [n09], [n10], [n11], [n12],
        [n13], [n14], [n15], [n16], [n17], [n18], [n19], [n20],
    ]
    write_submit_all(layout, ordered)
    write_pipeline_status(layout)

    print(f"Generated pipeline at: {layout.root}")
    print(f"  - submit_scripts/   ({sum(len(s) for s in ordered)} scripts)")
    print(f"  - submit_all.sh     (driver)")
    print(f"  - pipeline_status.sh (status check)")
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
