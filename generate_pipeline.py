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
  s02b_demux_bam     parallel demux that emits sorted BAMs (retains MM/ML
                     methylation tags from basecalling)
  s21_mijamp         MIJAMP preprocess per sample on the demux BAM + Bakta .fna

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


def barcode_filename_patterns(barcode: str) -> List[str]:
    """Return filename glob patterns to match a samplesheet barcode against
    dorado's actual output filenames. Handles both ONT kit naming (NB13) and
    dorado's internal naming (barcode13, zero-padded to two digits).

    For NB13 returns ['NB13', 'barcode13'].
    For barcode13 returns ['barcode13', 'NB13'].
    For NB5 returns ['NB5', 'barcode05'].
    """
    cands = [barcode]
    m = _re.match(r"^NB0*(\d+)$", barcode)
    if m:
        n = int(m.group(1))
        cands.append(f"barcode{n:02d}")
    m = _re.match(r"^barcode0*(\d+)$", barcode)
    if m:
        n = int(m.group(1))
        cands.append(f"NB{n}")
    # dedupe while preserving order
    seen = set()
    out = []
    for c in cands:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


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


def demux_classify_flag(cfg: dict) -> str:
    """Return the barcode-classification flag for the demux stages.

    Default is --no-classify, which trusts the BC tags dorado wrote during
    basecalling. With demux_reclassify: true, demux instead re-scans the
    reads against demux_kit_name (falling back to kit_name) — this is how you
    change the barcode kit without re-running the GPU-expensive basecall.
    """
    if cfg.get("demux_reclassify", False):
        kit = cfg.get("demux_kit_name") or cfg["kit_name"]
        return f"--kit-name {kit}"
    return "--no-classify"


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
            local ts line
            ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
            line="$(printf '%s [%s] %s %s' "$ts" "$level" "$JOB_NAME" "$*")"
            # Serialise the append. Dozens of s06 batch jobs write to this one
            # file concurrently, and appends are NOT atomic on NFS — without a
            # lock, lines get interleaved or silently dropped (which made
            # pipeline_status.sh report phantom "stalled" stages).
            if command -v flock >/dev/null 2>&1; then
                (
                    flock -w 15 9 2>/dev/null || true
                    printf '%s\\n' "$line" >> "$PIPELINE_LOG"
                ) 9>>"$PIPELINE_LOG.lock"
            else
                printf '%s\\n' "$line" >> "$PIPELINE_LOG"
            fi
            printf '%s\\n' "$line" >&2
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
                | grep -viE '(error[_ ]?rate|error[_ ]detection|error[_ ]adjustment|error[_ ]threshold|fraction[_ ]error|error[_ ]correction|nextDenovo|cns_align|Re-write workdir|empty job|log\\.critical|\\.tmpfrunt|sgs_fofn|hifi_fofn|Delete task:|plassembler|Unicycler has failed|no plasmids|uncorrected reads|Canu failed to correct|pymp-|_remove_temp_dir|multiprocessing/util|Device or resource busy|\\[(SCAN|WARN|OK|SKIP|START|DONE|FAIL|INFO|SUMMARY)\\])' \\
                | head -20 \\
                | sed -E 's/\\x1B\\[[0-9;]*[mGK]//g' \\
                | while IFS= read -r line; do
                    log_event SCAN "$line"
                done
        }}
        # samtools 0.1.x parses `sort -o out.bam -` as input=out.bam (in 0.1.x
        # -o is a boolean meaning "write to stdout"), so it tries to READ the
        # not-yet-existent output and dies with ENOENT, which in turn SIGPIPEs
        # minimap2 (rc=141). Catch that here rather than emitting confusing
        # "[bam_sort_core] fail to open file" errors.
        check_samtools() {{
            local envname="${{1:-this}}"
            if ! command -v samtools >/dev/null 2>&1; then
                log_event WARN "samtools not found on PATH in the '$envname' conda env"
                return 1
            fi
            local ver major
            ver=$(samtools --version 2>/dev/null | head -1 | awk '{{print $2}}')
            if [ -z "$ver" ]; then
                log_event WARN "samtools in '$envname' is pre-1.0 (no --version support). This pipeline needs samtools >= 1.10 — 0.1.x parses 'sort -o out.bam -' backwards. Fix with: mamba install -n $envname -c bioconda -c conda-forge 'samtools>=1.18'"
                return 1
            fi
            major=${{ver%%.*}}
            if ! [ "$major" -ge 1 ] 2>/dev/null; then
                log_event WARN "samtools $ver in '$envname' is too old; need >= 1.10. Fix with: mamba install -n $envname -c bioconda -c conda-forge 'samtools>=1.18'"
                return 1
            fi
            log_event OK "samtools $ver"
            return 0
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
    bams: Path
    mapping: Path
    meta_flye: Path
    bins: Path
    checkm2: Path
    gtdbtk: Path
    reports: Path
    mijamp: Path
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
            bams=root / "bams",
            mapping=root / "mapping",
            meta_flye=root / "meta_flye",
            bins=root / "bins",
            checkm2=root / "checkm2",
            gtdbtk=root / "gtdbtk",
            reports=root / "reports",
            mijamp=root / "mijamp",
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
            log_event OK "starting dorado demux ({demux_classify_flag(cfg)})"
            set +e
            {cfg['dorado_binary']} demux \\
                --output-dir {layout.fastq_raw} \\
                {demux_classify_flag(cfg)} \\
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
            # Dorado demux preserves the MinKNOW directory hierarchy
            # (experiment_id/sample_id/run_id/fastq_pass/barcodeNN/). Flatten
            # it so downstream stages can glob for *barcodeNN*.fastq directly
            # in {layout.fastq_raw}/.
            log_event OK "flattening nested demux output"
            find {layout.fastq_raw} -mindepth 2 -type f \\( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*summary*.tsv" \\) | while read f; do
                bn=$(basename "$f")
                if [ ! -e "{layout.fastq_raw}/$bn" ]; then
                    mv "$f" "{layout.fastq_raw}/$bn"
                fi
            done
            find {layout.fastq_raw} -mindepth 1 -type d -empty -delete 2>/dev/null || true
            log_event OK "demux complete"
        fi

        log_warn_errors_from "$STAGE_LOG"
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s02b demux to BAM (parallel branch off s01) ------------------------

def stage_demux_bam(cfg: dict, layout: Layout, hold: str) -> str:
    """A second demux that emits sorted BAMs (preserving the MM/ML
    methylation tags from basecalling) into <output>/bams/. Runs in parallel
    with s02_demux off the same s01_basecall output."""
    name = "s02b_demux_bam"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"],
        cfg.get("sge_h_vmem_demux_bam", "100G"),
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_demux"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += textwrap.dedent(
        f"""
        if ! check_samtools "{cfg['conda_env_demux']}"; then
            log_event WARN "aborting s02b — samtools >= 1.10 is required to sort/index the demuxed BAMs"
            exit 1
        fi

        mkdir -p {layout.bams}

        # Resume if any *.bam already exists in the demux BAM dir.
        if find {layout.bams} -type f -name "*.bam" 2>/dev/null | grep -q .; then
            log_event SKIP "demux BAM output already present in {layout.bams}"
        elif [ ! -s {layout.bam_path} ]; then
            log_event WARN "basecalls.bam missing, cannot demux to BAM"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
        else
            log_event OK "starting dorado demux (BAM, preserves MM/ML, {demux_classify_flag(cfg)})"
            set +e
            # NOTE: dorado's own --sort-bam is deliberately NOT used here. It
            # is incompatible with --kit-name (re-classification): the stage
            # exits 1 within seconds and prints nothing. Sorting with samtools
            # afterwards works in both modes and additionally gives us .bai
            # indexes, which dorado --sort-bam does not produce and which
            # modkit / MIJAMP need.
            {cfg['dorado_binary']} demux \\
                --output-dir {layout.bams} \\
                {demux_classify_flag(cfg)} \\
                --threads {cfg['sge_pe_parallel']} \\
                --emit-summary \\
                {layout.bam_path}
            rc=$?
            set -e
            if [ $rc -ne 0 ]; then
                log_event WARN "dorado demux (BAM) rc=$rc"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                exit $rc
            fi
            # Flatten the MinKNOW-style nested output. Dorado writes
            # bams/<experiment_id>/<sample_id>/<run_id>/bam_pass/barcodeNN/<file>.bam
            # — we move everything to bams/ so s21_mijamp can find it via glob.
            log_event OK "flattening nested BAM output"
            find {layout.bams} -mindepth 2 -type f \\( -name "*.bam" -o -name "*.bam.bai" -o -name "*summary*.tsv" \\) | while read f; do
                bn=$(basename "$f")
                if [ ! -e "{layout.bams}/$bn" ]; then
                    mv "$f" "{layout.bams}/$bn"
                fi
            done
            find {layout.bams} -mindepth 1 -type d -empty -delete 2>/dev/null || true

            # Sort + index each per-barcode BAM. Replaces dorado --sort-bam
            # (see note above) and produces the .bai files modkit/MIJAMP need.
            n_sorted=0
            n_sort_fail=0
            for b in {layout.bams}/*.bam; do
                [ -e "$b" ] || continue
                case "$b" in *.sorted.bam) continue ;; esac
                if [ -s "${{b}}.bai" ]; then
                    continue   # already sorted+indexed by an earlier run
                fi
                set +e
                samtools sort -@ {cfg['sge_pe_parallel']} -o "${{b}}.sorting" "$b" \\
                    && mv -f "${{b}}.sorting" "$b" \\
                    && samtools index "$b"
                src=$?
                set -e
                if [ $src -ne 0 ]; then
                    log_event WARN "samtools sort/index failed for $(basename "$b") rc=$src"
                    rm -f "${{b}}.sorting"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    n_sort_fail=$((n_sort_fail+1))
                else
                    n_sorted=$((n_sorted+1))
                fi
            done

            n_bams=$(find {layout.bams} -maxdepth 1 -type f -name "*.bam" | wc -l)
            log_event OK "demux (BAM) produced $n_bams files (sorted+indexed: $n_sorted, failed: $n_sort_fail)"
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
        # Build a multi-pattern find that accepts both NB13 and barcode13
        # naming so the samplesheet can use either convention.
        pats = barcode_filename_patterns(s.barcode)
        name_clauses = " -o ".join(
            f'-name "*{p}*.fastq" -o -name "*{p}*.fastq.gz"' for p in pats
        )
        body += textwrap.dedent(
            f"""
            # {s.sample} <- {s.barcode} (matches: {", ".join(pats)})
            if [ -s {out} ]; then
                log_event SKIP "{s.sample} (fastq exists)"
            else
                files=$(find {layout.fastq_raw} -type f \\( {name_clauses} \\) 2>/dev/null)
                if [ -z "$files" ]; then
                    log_event WARN "{s.sample}: no demux output matching any of: {', '.join(pats)}"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    : > {out}
                    for f in $files; do
                        case "$f" in
                            *.gz) zcat "$f" >> {out} ;;
                            *)    cat  "$f" >> {out} ;;
                        esac
                    done
                    # Report reads/bases/depth now rather than letting a
                    # low-coverage sample fail opaquely at s05_subsample.
                    nreads=$(awk 'NR%4==1' {out} | wc -l)
                    nbases=$(awk 'NR%4==2 {{s+=length($0)}} END {{print s+0}}' {out})
                    depth=$(awk -v b="$nbases" -v g={s.genome_length} \\
                        'BEGIN {{printf "%.2f", (g>0 ? b/g : 0)}}')
                    log_event OK "{s.sample} reads=$nreads bases=$nbases depth=${{depth}}x"
                    if awk -v d="$depth" 'BEGIN {{exit !(d < 25)}}'; then
                        log_event WARN "{s.sample}: depth ${{depth}}x is below autocycler's 25x minimum — this sample WILL fail at s05_subsample. Check that barcode {s.barcode} actually carried library material."
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    fi
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


# ---- s18 annotate (Bakta for Bacteria, Prokka for Archaea) --------------

def stage_annotate(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    """Annotate every assembly, choosing BOTH the annotator and the --proteins
    reference from GTDB-Tk's classification.

        Bacteria -> bakta  --proteins <RefSeq type-strain proteome>
        Archaea  -> prokka --kingdom Archaea --proteins <same>

    The per-genome decisions are made by _resolve_annotation_plan.py, which
    runs first inside this stage and writes annotation_plan.tsv. The
    type-strain proteomes come from a LOCAL directory built offline by
    build_type_strain_db.py (the cluster has no internet access).
    """
    name = "s18_annotate"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    # Bakta and Prokka live in the same conda env.
    body += conda_block(cfg["conda_env_bakta"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"

    # Standalone resolver, written once next to the submit scripts.
    plan_py = layout.submit_dir / "_resolve_annotation_plan.py"
    plan_py.write_text(ANNOTATION_PLAN_PY)
    plan_py.chmod(plan_py.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    # Curated samplesheet metadata, for the isolate consensuses.
    meta_path = layout.annotate / "sample_meta.tsv"
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    with meta_path.open("w") as fh:
        fh.write("sample\tgenus\tspecies\tstrain\tkingdom\tgram\tplasmid\treference_proteins\n")
        for s in samples:
            fh.write("\t".join([
                s.sample, s.genus, s.species, s.strain, s.kingdom,
                s.gram, s.plasmid, s.reference_proteins,
            ]) + "\n")

    plan     = f"{layout.annotate}/annotation_plan.tsv"
    missing  = f"{layout.annotate}/missing_species.txt"
    resolver_out = f"{layout.annotate}/plan_resolver.out"
    min_len  = int(cfg.get("bakta_min_contig_length", 200))
    threads  = int(cfg.get("bakta_threads", 20))
    ts_db    = cfg.get("type_strain_db_dir", "") or ""

    # Candidate locations searched on the cluster at runtime when
    # type_strain_db_dir is unset (or points somewhere unusable). The first
    # one holding a valid index.tsv is adopted, so a DB built and rsync'd
    # earlier gets picked up without editing the config.
    search_paths: List[str] = [f"{layout.root}/type_strain_db"]
    search_paths += [str(p) for p in cfg.get("type_strain_db_search_paths", []) or []]
    ts_search = os.pathsep.join(search_paths)

    body += textwrap.dedent(
        f"""
        mkdir -p {layout.annotate}

        log_event OK "resolving annotation plan from GTDB-Tk classifications"
        set +e
        python3 {plan_py} \\
            --gtdbtk-dir {layout.gtdbtk}/output \\
            --consensus-dir {layout.consensus} \\
            --bins-dir {layout.bins} \\
            --annotate-dir {layout.annotate} \\
            --type-strain-db "{ts_db}" \\
            --type-strain-search "{ts_search}" \\
            --sample-meta {meta_path} \\
            --out-plan {plan} \\
            --out-missing {missing} \\
            > {resolver_out} 2>&1
        rc=$?
        set -e
        cat {resolver_out} >&2
        if [ $rc -ne 0 ] || [ ! -s {plan} ]; then
            log_event WARN "could not build annotation plan (rc=$rc)"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            exit 1
        fi
        # Surface which type-strain DB was actually used (config vs discovered).
        ts_line=$(grep -m1 '^type_strain_db:' {resolver_out} || true)
        [ -n "$ts_line" ] && log_event OK "$ts_line"
        if grep -q '^type_strain_db: none found' {resolver_out}; then
            log_event WARN "no type-strain DB found — annotating WITHOUT --proteins. Build one with build_type_strain_db.py (see README) and drop it at {layout.root}/type_strain_db, or set type_strain_db_dir."
        fi
        log_event OK "annotation plan: $(tail -n +2 {plan} | wc -l) genomes"
        if [ -s {missing} ]; then
            log_event WARN "$(wc -l < {missing}) species have no local type-strain proteome. See {missing} — run build_type_strain_db.py on an internet-connected machine, rsync the result over, and re-run this stage."
        fi

        # Walk the plan. Process substitution (not a pipe) so the loop runs in
        # this shell and STAGE_FAIL_COUNT survives.
        while IFS=$'\\t' read -r gid fasta atype domain tool genus epithet strain gram plasmid gtdb_sp proteins psrc out_dir; do
            [ -n "$gid" ] || continue
            # The plan writes __NONE__ for empty fields so that consecutive
            # tabs never collapse (bash treats tab as IFS-whitespace and would
            # otherwise shift every subsequent column left). Undo that here.
            for _v in gid fasta atype domain tool genus epithet strain gram \\
                      plasmid gtdb_sp proteins psrc out_dir; do
                # Belt and braces: strip any stray carriage return (a CRLF
                # plan would otherwise leave \\r on the final column and we
                # would silently create paths ending in a control character).
                printf -v "$_v" '%s' "${{!_v%$'\\r'}}"
                [ "${{!_v}}" = "__NONE__" ] && printf -v "$_v" '%s' ""
            done

            if [ "$tool" = "prokka" ]; then
                marker="$out_dir/${{gid}}.gff"
            else
                marker="$out_dir/${{gid}}.gff3"
            fi

            if [ -s "$marker" ]; then
                log_event SKIP "$gid ($tool output exists)"
                continue
            fi
            if [ ! -s "$fasta" ]; then
                log_event WARN "$gid: input fasta missing ($fasta)"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                continue
            fi

            # --proteins only if the resolved file actually exists on disk
            prot_flag=""
            if [ -n "$proteins" ] && [ -s "$proteins" ]; then
                prot_flag="--proteins $proteins"
            elif [ -n "$proteins" ]; then
                log_event WARN "$gid: proteins file listed but missing ($proteins), annotating without it"
            fi

            genus_flag=""
            [ -n "$genus" ]   && genus_flag="--genus $genus"
            species_flag=""
            [ -n "$epithet" ] && species_flag="--species $epithet"
            strain_flag=""
            [ -n "$strain" ]  && strain_flag="--strain $strain"

            rm -rf "$out_dir"
            mkdir -p "$out_dir"
            set +e
            if [ "$tool" = "prokka" ]; then
                log_event OK "$gid: prokka (Archaea, ${{gtdb_sp:-unclassified}}, proteins=$psrc)"
                prokka \\
                    --outdir "$out_dir" \\
                    --prefix "$gid" \\
                    --locustag "$gid" \\
                    --kingdom Archaea \\
                    $genus_flag $species_flag $strain_flag \\
                    $prot_flag \\
                    --cpus {threads} \\
                    --force \\
                    "$fasta"
                arc=$?
            else
                gram_flag=""
                case "$gram" in
                    pos) gram_flag="--gram +" ;;
                    neg) gram_flag="--gram -" ;;
                    unknown) gram_flag="--gram ?" ;;
                esac
                plasmid_flag=""
                [ -n "$plasmid" ] && plasmid_flag="--plasmid $plasmid"
                log_event OK "$gid: bakta (${{domain:-unknown}}, ${{gtdb_sp:-unclassified}}, proteins=$psrc)"
                bakta \\
                    --db {cfg['bakta_db']} \\
                    --min-contig-length {min_len} \\
                    --prefix "$gid" \\
                    --output "$out_dir" \\
                    $genus_flag $species_flag $strain_flag \\
                    $gram_flag $plasmid_flag \\
                    $prot_flag \\
                    --threads {threads} \\
                    --force \\
                    "$fasta"
                arc=$?
            fi
            set -e

            # A malformed reference proteome must never cost us the whole
            # annotation. If the run failed while --proteins was in play,
            # retry once without it — an annotation without the extra
            # evidence is far better than none.
            if {{ [ $arc -ne 0 ] || [ ! -s "$marker" ]; }} && [ -n "$prot_flag" ]; then
                log_event WARN "$gid: $tool failed (rc=$arc) with --proteins; retrying without it"
                rm -rf "$out_dir"; mkdir -p "$out_dir"
                set +e
                if [ "$tool" = "prokka" ]; then
                    prokka \\
                        --outdir "$out_dir" \\
                        --prefix "$gid" \\
                        --locustag "$gid" \\
                        --kingdom Archaea \\
                        $genus_flag $species_flag $strain_flag \\
                        --cpus {threads} \\
                        --force \\
                        "$fasta"
                    arc=$?
                else
                    bakta \\
                        --db {cfg['bakta_db']} \\
                        --min-contig-length {min_len} \\
                        --prefix "$gid" \\
                        --output "$out_dir" \\
                        $genus_flag $species_flag $strain_flag \\
                        $gram_flag $plasmid_flag \\
                        --threads {threads} \\
                        --force \\
                        "$fasta"
                    arc=$?
                fi
                set -e
                if [ $arc -eq 0 ] && [ -s "$marker" ]; then
                    log_event WARN "$gid annotated by $tool WITHOUT --proteins (the reference at $proteins was rejected — check its format)"
                fi
            fi

            if [ $arc -ne 0 ] || [ ! -s "$marker" ]; then
                log_event WARN "$gid: $tool failed (rc=$arc)"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                log_event OK "$gid annotated by $tool"
            fi
        done < <(tail -n +2 {plan})

        log_warn_errors_from "$STAGE_LOG"
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s13 map reads (minimap2 + samtools) ---------------------------------

def stage_map_reads(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s13_map_reads"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"],
        cfg.get("sge_h_vmem_mapping", "100G"),
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_mapping"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    threads = int(cfg.get("minimap2_threads", 20))
    body += textwrap.dedent(
        f"""
        if ! check_samtools "{cfg['conda_env_mapping']}"; then
            log_event WARN "aborting s13 — a working samtools >= 1.10 is required"
            exit 1
        fi
        if ! command -v minimap2 >/dev/null 2>&1; then
            log_event WARN "minimap2 not found in the '{cfg['conda_env_mapping']}' env"
            exit 1
        fi

        """
    )
    body += f"mkdir -p {layout.mapping}\n\n"
    for s in samples:
        sdir = f"{layout.mapping}/{s.sample}"
        marker = f"{sdir}/flagstat.txt"
        # Resume marker is a .done sentinel, NOT the presence of unmapped.fastq:
        # a deeply-sequenced isolate can legitimately have ZERO unmapped reads,
        # which would make an "-s unmapped.fastq" guard re-run forever.
        done_marker = f"{sdir}/.mapping_done"
        body += textwrap.dedent(
            f"""
            if ! mkdir -p {sdir} 2>/dev/null || [ ! -d {sdir} ]; then
                log_event WARN "{s.sample}: cannot create {sdir} (disk full / quota / permissions?)"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            elif [ -f {done_marker} ] && [ -s {marker} ]; then
                log_event SKIP "{s.sample} (already mapped)"
            elif [ ! -s {layout.consensus}/{s.sample}.fasta ]; then
                log_event WARN "{s.sample}: consensus fasta missing, skipping"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            elif [ ! -s {layout.fastq_filtered}/{s.sample}.fastq ]; then
                log_event WARN "{s.sample}: filtered reads missing, skipping"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                set +e
                # pipefail so a minimap2 failure is not masked by samtools sort
                # exiting 0. minimap2 stderr is kept (previously /dev/null'd,
                # which hid the real cause of failures).
                set -o pipefail
                minimap2 -ax map-ont -t {threads} \\
                    {layout.consensus}/{s.sample}.fasta \\
                    {layout.fastq_filtered}/{s.sample}.fastq \\
                    2> {sdir}/minimap2.stderr \\
                    | samtools sort -@ {threads} -o {sdir}/mapped.bam -
                rc=$?
                set +o pipefail
                # Never trust rc alone — verify the artifact actually exists.
                # samtools has been observed exiting 0 after failing to create
                # its output file, which previously produced a bogus [OK].
                if [ $rc -ne 0 ] || [ ! -s {sdir}/mapped.bam ]; then
                    log_event WARN "{s.sample}: mapping failed (rc=$rc, mapped.bam missing/empty)"
                    if [ -s {sdir}/minimap2.stderr ]; then
                        log_event WARN "{s.sample}: minimap2 said: $(tail -3 {sdir}/minimap2.stderr | tr '\\n' ' ')"
                    fi
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    samtools index {sdir}/mapped.bam
                    samtools flagstat {sdir}/mapped.bam > {marker}
                    # -f 4 selects unmapped records
                    samtools fastq -f 4 {sdir}/mapped.bam > {sdir}/unmapped.fastq 2>/dev/null
                    mapped=$(grep -E '[0-9]+ \\+ [0-9]+ mapped \\(' {marker} | head -1 | awk '{{print $1}}')
                    total=$(head -1 {marker} | awk '{{print $1}}')
                    unmapped_reads=$(awk 'NR%4==1' {sdir}/unmapped.fastq 2>/dev/null | wc -l)
                    # Guard against an unparseable/empty flagstat
                    if [ -z "$mapped" ] || [ -z "$total" ]; then
                        log_event WARN "{s.sample}: could not parse flagstat ({marker} empty or malformed)"
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    else
                        pct=$(awk -v m="$mapped" -v t="$total" \\
                            'BEGIN {{printf "%.2f", (t>0 ? 100*m/t : 0)}}')
                        log_event OK "{s.sample} mapped=$mapped/$total (${{pct}}%) unmapped_reads=$unmapped_reads"
                        touch {done_marker}
                    fi
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
    body += textwrap.dedent(
        f"""
        if ! check_samtools "{cfg['conda_env_metabat']}"; then
            log_event WARN "aborting s15 — samtools >= 1.10 is required to build the depth BAM"
            exit 1
        fi

        """
    )
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
    ck_tmp = cfg.get("checkm2_tmpdir", "") or ""
    body += textwrap.dedent(
        f"""

        # Keep python multiprocessing's temp dir OFF NFS. Deleting a
        # still-open file on NFS becomes a "silly rename" that cannot be
        # unlinked, so CheckM2's shutdown finaliser throws
        #   OSError: [Errno 16] Device or resource busy: .../pymp-XXXX
        # after it has already written quality_report.tsv. Results are fine,
        # but it is noisy and leaves stale pymp-* dirs on scratch.
        CK_TMP="{ck_tmp}"
        [ -n "$CK_TMP" ] || CK_TMP="/tmp/checkm2.$$"
        if mkdir -p "$CK_TMP" 2>/dev/null; then
            export TMPDIR="$CK_TMP"
            log_event OK "TMPDIR=$TMPDIR (node-local, avoids NFS unlink errors)"
        else
            log_event WARN "could not create $CK_TMP; leaving TMPDIR=${{TMPDIR:-/tmp}} (expect benign pymp-* unlink tracebacks)"
            CK_TMP=""
        fi

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

        # Tidy our node-local temp dir (and any pymp-* left behind on it).
        if [ -n "${{CK_TMP:-}}" ] && [ -d "$CK_TMP" ]; then
            rm -rf "$CK_TMP" 2>/dev/null || true
        fi

        log_warn_errors_from "$STAGE_LOG"
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s18 GTDB-Tk classification (consensuses + bins) --------------------

def stage_classify(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    name = "s16_classify"
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

        # Build the GTDB-Tk batchfile from the Autocycler consensuses and the
        # MetaBAT2 bins DIRECTLY — not from Bakta's .fna.
        #
        # GTDB-Tk calls genes itself (Prodigal) and only needs nucleotide
        # FASTA, so depending on Bakta bought nothing while forcing GTDB-Tk to
        # queue behind ~1h of annotation and making classification fail
        # whenever annotation did.
        #
        # Format: <fasta_path>\\t<genome_id>
        # Bin IDs are <sample>_bin.N so origin sample is traceable in the
        # GTDB-Tk summary rows.
        : > {batchfile}
        """
    )
    for s in samples:
        cons = f"{layout.consensus}/{s.sample}.fasta"
        body += textwrap.dedent(
            f"""
            if [ -s {cons} ]; then
                printf '%s\\t%s\\n' '{cons}' '{s.sample}' >> {batchfile}
            else
                log_event WARN "{s.sample}: missing consensus ({cons}), excluded from GTDB-Tk"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            fi
            """
        )
    body += textwrap.dedent(
        f"""

        # Add every MetaBAT2 bin (raw output — no annotation needed).
        for f in {layout.bins}/*/*.fa; do
            [ -e "$f" ] || continue
            genome_id=$(basename "$f" .fa)
            printf '%s\\t%s\\n' "$f" "$genome_id" >> {batchfile}
        done

        n_entries=$(wc -l < {batchfile})
        log_event OK "GTDB-Tk batchfile has $n_entries entries (consensuses + bins)"
        if [ "$n_entries" -eq 0 ]; then
            log_event WARN "no assemblies found, skipping GTDB-Tk entirely"
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
            \"\"\"Feature counts from whichever annotator ran.

            Bakta writes <prefix>.json (preferred) and <prefix>.txt.
            Prokka (used for Archaea) writes only <prefix>.txt, in the form
                CDS: 4200
                tRNA: 86
            Both live at <annotate_dir>/<genome_id>/<genome_id>.*
            \"\"\"
            sample_dir = p.parent
            stem = p.stem
            counts = {{"CDS": 0, "tRNA": 0, "rRNA": 0, "tmRNA": 0,
                      "ncRNA": 0, "CRISPR": 0, "sORF": 0, "oriC": 0, "oriT": 0}}
            j = sample_dir / f"{{stem}}.json"
            if j.exists():
                d = load_bakta_json(j)
                for f in d.get("features", []):
                    t = f.get("type", "")
                    if t in counts:
                        counts[t] += 1
                    elif t == "ncRNA-region":
                        counts["ncRNA"] += 1
                if any(counts.values()):
                    return counts
            # Prokka (or Bakta without JSON): parse the .txt summary.
            t = sample_dir / f"{{stem}}.txt"
            if t.exists():
                # Prokka keys: CDS, rRNA, tRNA, tmRNA, misc_RNA, repeat_region
                alias = {{"misc_RNA": "ncRNA", "repeat_region": "CRISPR"}}
                for line in t.read_text(errors="ignore").splitlines():
                    if ":" not in line:
                        continue
                    k, _, v = line.partition(":")
                    k = k.strip()
                    k = alias.get(k, k)
                    if k in counts:
                        try:
                            counts[k] = int(v.strip().split()[0])
                        except Exception:
                            pass
            return counts

        def annotator_of(genome_id):
            \"\"\"Which tool annotated this genome, from annotation_plan.tsv.\"\"\"
            plan = ROOT / "bakta" / "annotation_plan.tsv"
            if not plan.exists():
                return ""
            try:
                with plan.open() as fh:
                    for row in csv.DictReader(fh, delimiter="\\t"):
                        if (row.get("genome_id") or "").strip() == genome_id:
                            return (row.get("tool") or "").strip()
            except Exception:
                pass
            return ""

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

        def load_plassembler_clusters():
            \"\"\"Parse stdout/s09_cluster.log to count qc_pass clusters that
            had plassembler-derived source contigs, per sample. Used to decide
            whether a consensus should be split into chromosome + plasmid(s)
            in the report. Without Plassembler validation, multi-contig
            consensuses are treated as one main assembly.\"\"\"
            log_path = ROOT / "stdout" / "s09_cluster.log"
            result = {{}}
            if not log_path.exists():
                return result
            current_sample = None
            in_cluster = False
            has_plass = False
            for line in log_path.read_text(errors="ignore").splitlines():
                stripped = line.strip()
                m = re.search(r'--autocycler_dir\\s+\\S+/([^/\\s]+)\\s*$', stripped)
                if m:
                    current_sample = m.group(1)
                    result.setdefault(current_sample, 0)
                    in_cluster = False
                    has_plass = False
                    continue
                if stripped.startswith("Cluster "):
                    in_cluster = True
                    has_plass = False
                    continue
                if not in_cluster:
                    continue
                if "_plassembler_" in stripped:
                    has_plass = True
                elif stripped.startswith("passed QC"):
                    if has_plass and current_sample:
                        result[current_sample] = result.get(current_sample, 0) + 1
                    in_cluster = False
                    has_plass = False
                elif stripped.startswith("failed QC"):
                    in_cluster = False
                    has_plass = False
            return result

        checkm2     = load_checkm2()
        gtdbtk      = load_gtdbtk()
        expected    = load_expected()
        plass_clust = load_plassembler_clusters()

        rows = []
        for s in SAMPLES:
            raw_reads  = fastq_count(ROOT / "fastq" / f"{{s}}.fastq")
            filt_reads = fastq_count(ROOT / "fastq_filtered" / f"{{s}}.fastq")
            flag = parse_flagstat(ROOT / "mapping" / s / "flagstat.txt")
            total, mapped = (flag if flag else (0, 0))
            unmapped_reads = fastq_count(ROOT / "mapping" / s / "unmapped.fastq")
            cons = ROOT / "consensus" / f"{{s}}.fasta"
            c_n, c_len, c_long = fasta_stats(cons)
            # Split consensus into chromosome (largest contig) + plasmid contigs,
            # but ONLY treat non-largest contigs as plasmids if Plassembler
            # validated at least one plasmid cluster for this sample (parsed
            # from the s09_cluster log). Otherwise the consensus is shown as a
            # single main assembly even if it has multiple contigs.
            contigs_sorted = fasta_per_contig(cons)
            has_plassembler_plasmid = plass_clust.get(s, 0) > 0
            if contigs_sorted:
                chrom_name, chrom_len = contigs_sorted[0]
            else:
                chrom_name, chrom_len = "", 0
            if has_plassembler_plasmid and len(contigs_sorted) > 1:
                plasmid_contigs = contigs_sorted[1:]
            else:
                plasmid_contigs = []
            plasmid_count = len(plasmid_contigs)
            plasmid_total_bp = sum(L for _, L in plasmid_contigs)
            plasmid_names = ";".join(n for n, _ in plasmid_contigs)
            plasmid_lengths = ";".join(str(L) for _, L in plasmid_contigs)
            n_plassembler_clusters = plass_clust.get(s, 0)
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
                "plassembler_clusters": n_plassembler_clusters,
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
                "annotator": annotator_of(s),
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
                # Bins are annotated too (Bakta or Prokka depending on domain)
                bfeats = bakta_feature_counts(
                    ROOT / "bakta" / bid / f"{{bid}}.fna")
                rows.append({{
                    "sample": s,
                    "assembly_type": "metagenomic_bin",
                    "genome_id": bid,
                    "annotator": annotator_of(bid),
                    "cds": bfeats.get("CDS", 0),
                    "tRNA": bfeats.get("tRNA", 0),
                    "rRNA": bfeats.get("rRNA", 0),
                    "tmRNA": bfeats.get("tmRNA", 0),
                    "ncRNA": bfeats.get("ncRNA", 0),
                    "CRISPR": bfeats.get("CRISPR", 0),
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
                  "plasmid_lengths","plassembler_clusters",
                  "meta_assembly_contigs",
                  "meta_assembly_length_bp","n_bins","completeness",
                  "contamination","classification","domain","phylum","class",
                  "order","family","genus","species","expected_genus",
                  "expected_species","genus_match","species_match",
                  "annotator","cds","tRNA","rRNA","tmRNA","ncRNA","CRISPR"]
        with OUT.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\\t",
                               extrasaction="ignore", lineterminator="\\n")
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
ANNOTATION_PLAN_PY = r'''#!/usr/bin/env python3
"""Decide, per genome, which annotator to run and which --proteins to pass.

Runs on the cluster AFTER GTDB-Tk. For every assembly (Autocycler consensus
or MetaBAT2 bin):
  * read the GTDB-Tk classification
  * Bacteria -> bakta,  Archaea -> prokka
  * resolve the species to a local RefSeq type-strain proteome (offline)

Writes annotation_plan.tsv for the annotation stage to consume, plus
missing_species.txt listing anything that could not be resolved (feed that to
build_type_strain_db.py on an internet-connected machine).
"""
import argparse, csv, os, re, sys
from pathlib import Path


def pick_type_strain_db(explicit, search_paths):
    """Choose which type-strain DB to use, at RUNTIME on the cluster.

    An explicitly configured directory wins if it is usable. Otherwise the
    first candidate that contains a non-empty index.tsv is adopted, so a DB
    that was built and rsync'd earlier is picked up automatically without
    anyone having to edit the config.

    Returns (path, how) where how is one of: config, auto, none.
    """
    def usable(p):
        if not p:
            return None
        q = Path(os.path.expandvars(os.path.expanduser(str(p))))
        idx = q / "index.tsv"
        try:
            if idx.is_file() and idx.stat().st_size > 0:
                return str(q)
        except OSError:
            pass
        return None

    if explicit:
        got = usable(explicit)
        if got:
            return got, "config"
        sys.stderr.write(
            f"warning: type_strain_db_dir={explicit!r} has no usable "
            f"index.tsv; falling back to auto-discovery\n")

    for cand in search_paths:
        got = usable(cand)
        if got:
            return got, "auto"
    return "", "none"


def normalize_taxon(name):
    """Strip GTDB placeholder suffixes (Methanobrevibacter_A -> Methanobrevibacter)."""
    if not name:
        return ""
    name = re.sub(r"^[a-z]__", "", name.strip())
    parts = [re.sub(r"_[A-Z]+$", "", p) for p in name.split()]
    return " ".join(p for p in parts if p).strip()


def species_key(species):
    return " ".join(normalize_taxon(species).lower().split())


def parse_classification(cls):
    out = {}
    for field in (cls or "").split(";"):
        field = field.strip()
        for pref, name in (("d__", "domain"), ("p__", "phylum"), ("c__", "class"),
                           ("o__", "order"), ("f__", "family"), ("g__", "genus"),
                           ("s__", "species")):
            if field.startswith(pref):
                out[name] = field[len(pref):].strip()
                break
    return out


def load_gtdbtk(gtdbtk_dir):
    """genome_id -> {domain, genus, species, classification}"""
    res = {}
    d = Path(gtdbtk_dir)
    for f in sorted(d.glob("gtdbtk.*.summary.tsv")):
        # bac120 -> Bacteria, ar53 -> Archaea (also read from classification)
        with f.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gid = (row.get("user_genome") or "").strip()
                if not gid:
                    continue
                cls = row.get("classification", "") or ""
                t = parse_classification(cls)
                domain = t.get("domain", "")
                if not domain:
                    domain = "Archaea" if "ar53" in f.name else "Bacteria"
                res[gid] = {
                    "domain": domain,
                    "genus": t.get("genus", ""),
                    "species": t.get("species", ""),
                    "classification": cls,
                }
    return res


def load_type_strain_index(db_dir):
    """species_key -> abs proteome path;  plus genus -> abs proteome path"""
    by_species, by_genus = {}, {}
    if not db_dir:
        return by_species, by_genus
    d = Path(db_dir)
    idx = d / "index.tsv"
    if not idx.exists():
        sys.stderr.write(f"warning: no index.tsv in {db_dir}\n")
        return by_species, by_genus
    with idx.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rel = (row.get("proteome_rel_path") or "").strip()
            if not rel:
                continue
            p = d / rel
            if not p.exists() or p.stat().st_size == 0:
                continue
            sk = (row.get("species_key") or "").strip().lower()
            gn = (row.get("genus") or "").strip().lower()
            if sk:
                by_species[sk] = str(p)
            if gn and gn not in by_genus:
                by_genus[gn] = str(p)
    return by_species, by_genus


def resolve_proteins(species, genus, by_species, by_genus):
    """Cascade: exact species -> normalised species -> genus -> none."""
    if species:
        raw = " ".join(species.lower().split())
        if raw in by_species:
            return by_species[raw], "species"
        norm = species_key(species)
        if norm and norm in by_species:
            return by_species[norm], "species-normalised"
    if genus:
        g = normalize_taxon(genus).lower()
        if g in by_genus:
            return by_genus[g], "genus-fallback"
    return "", "none"


def load_sample_meta(path):
    """sample -> dict of curated samplesheet metadata."""
    meta = {}
    if not path:
        return meta
    p = Path(path)
    if not p.exists():
        return meta
    with p.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s = (row.get("sample") or "").strip()
            if s:
                meta[s] = {k: (v or "").strip() for k, v in row.items()}
    return meta


def epithet_of(species):
    """'Escherichia coli' -> 'coli'.  Bakta/Prokka want genus and epithet
    as separate flags."""
    toks = normalize_taxon(species).split()
    return " ".join(toks[1:]) if len(toks) >= 2 else ""


GENBANK_SUFFIXES = {".gb", ".gbk", ".gbf", ".gbff", ".genbank", ".embl", ".dat"}


def prepare_proteins(src, tool, workdir):
    """Return a --proteins path that the chosen annotator will actually accept.

    Bakta's parser (bakta/expert/protein_sequences.py) is:

        cols = record.description.split(' ', 1)[1].split('~~~')

    i.e. it splits the header at the FIRST SPACE and then splits the
    remainder on '~~~', expecting three fields:

        >SeqID<space>gene~~~product~~~dbxrefs

    Anything else is discarded per-record ('wrong description format in
    FASTA user protein file'), and a header containing no space at all
    raises IndexError and aborts the whole run:

        ERROR - EXPERT_AA_SEQ - provided user proteins file Fasta format not valid!
        IndexError: list index out of range

    A plain RefSeq header ('>WP_000123456.1 chorismate mutase
    [Escherichia coli]') therefore has to be rewritten. GenBank input is
    parsed natively by Bakta — and by Prokka — so it passes straight
    through and is the format to prefer; Prokka also tolerates plain FASTA.

    Converted files are cached in workdir, so this costs one pass the first
    time and nothing on resume.
    """
    if not src:
        return "", "none"
    p = Path(src)
    try:
        if not p.is_file() or p.stat().st_size == 0:
            return "", "missing"
    except OSError:
        return "", "missing"

    # GenBank/EMBL: both tools read it directly.
    if p.suffix.lower() in GENBANK_SUFFIXES:
        return str(p), "genbank"
    # Prokka copes with plain FASTA (description becomes /product).
    if tool != "bakta":
        return str(p), "asis"

    first = ""
    with p.open(errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                first = line[1:].rstrip("\n")
                break
    if not first:
        return "", "empty"
    # Already in Bakta's shape? Needs a space, then exactly 3 '~~~' fields.
    head, sep, rest = first.partition(" ")
    if sep and len(rest.split("~~~")) == 3:
        return str(p), "asis"

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / (p.stem + ".bakta.faa")
    if out.is_file() and out.stat().st_size > 0:
        return str(out), "converted-cached"

    n = 0
    tmp = out.with_suffix(".tmp")
    with p.open(errors="ignore") as fin, tmp.open("w") as fo:
        for line in fin:
            if line.startswith(">"):
                hdr = line[1:].strip()
                sid, _, desc = hdr.partition(" ")
                # drop a trailing '[Organism name]' and RefSeq's
                # 'MULTISPECIES: ' prefix, neither of which is a product name
                desc = re.sub(r"\s*\[[^\]]*\]\s*$", "", desc)
                desc = re.sub(r"^MULTISPECIES:\s*", "", desc).strip()
                # '~~~' inside a product would corrupt the field split
                desc = desc.replace("~~~", " ")
                if not desc:
                    desc = "hypothetical protein"
                # >SeqID<space>gene~~~product~~~dbxrefs   (gene left empty)
                fo.write(f">{sid} ~~~{desc}~~~RefSeq:{sid}\n")
                n += 1
            else:
                fo.write(line)
    if n == 0:
        tmp.unlink(missing_ok=True)
        return "", "empty"
    tmp.replace(out)
    return str(out), f"converted({n})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtdbtk-dir", required=True)
    ap.add_argument("--consensus-dir", required=True)
    ap.add_argument("--bins-dir", required=True)
    ap.add_argument("--annotate-dir", required=True)
    ap.add_argument("--type-strain-db", default="",
                    help="explicit type-strain DB dir (from pipeline_config.yaml)")
    ap.add_argument("--type-strain-search", default="",
                    help="os.pathsep-separated candidate dirs to auto-discover "
                         "an existing DB when --type-strain-db is unset/unusable")
    ap.add_argument("--sample-meta", default="",
                    help="TSV of curated samplesheet metadata for consensuses")
    ap.add_argument("--out-plan", required=True)
    ap.add_argument("--out-missing", required=True)
    args = ap.parse_args()

    search = [p for p in args.type_strain_search.split(os.pathsep) if p.strip()]
    db_dir, how = pick_type_strain_db(args.type_strain_db, search)
    if how == "config":
        print(f"type_strain_db: {db_dir} (from config)")
    elif how == "auto":
        print(f"type_strain_db: {db_dir} (auto-discovered existing DB)")
    else:
        print("type_strain_db: none found — annotating without --proteins")

    tax = load_gtdbtk(args.gtdbtk_dir)
    by_species, by_genus = load_type_strain_index(db_dir)
    meta = load_sample_meta(args.sample_meta)
    sys.stderr.write(
        f"type-strain index: {len(by_species)} species, {len(by_genus)} genera\n")

    genomes = []
    for f in sorted(Path(args.consensus_dir).glob("*.fasta")):
        genomes.append((f.stem, str(f), "isolate_consensus"))
    bins_root = Path(args.bins_dir)
    if bins_root.exists():
        for f in sorted(bins_root.glob("*/*.fa")):
            genomes.append((f.stem, str(f), "metagenomic_bin"))

    rows, missing = [], []
    for gid, path, kind in genomes:
        t = tax.get(gid, {})
        domain = t.get("domain", "")
        gtdb_species = t.get("species", "")
        gtdb_genus = t.get("genus", "")
        m = meta.get(gid, {})

        # Archaea -> Prokka, everything else -> Bakta. An unclassified genome
        # falls through to Bakta, the safer general-purpose default.
        tool = "prokka" if domain.lower().startswith("archae") else "bakta"

        # Organism metadata: curated samplesheet values win for the isolate
        # consensuses (they are deliberate); bins have no samplesheet entry so
        # they take GTDB-Tk's assignment.
        if kind == "isolate_consensus" and (m.get("genus") or m.get("species")):
            genus = m.get("genus", "")
            epithet = m.get("species", "")
            strain = m.get("strain", "") or gid
        else:
            genus = normalize_taxon(gtdb_genus)
            epithet = epithet_of(gtdb_species)
            strain = gid

        # --proteins: an explicit samplesheet reference_proteins is a curator
        # override and wins; otherwise use the RefSeq type strain for whatever
        # species GTDB-Tk assigned.
        override = m.get("reference_proteins", "")
        if override:
            proteins, src = override, "samplesheet"
        else:
            proteins, src = resolve_proteins(
                gtdb_species, gtdb_genus, by_species, by_genus)
            if not proteins and (gtdb_species or gtdb_genus):
                missing.append(gtdb_species or gtdb_genus)

        # Make the file digestible by whichever annotator will read it.
        if proteins:
            prepared, how = prepare_proteins(
                proteins, tool, Path(args.annotate_dir) / "proteins_prepared")
            if not prepared:
                sys.stderr.write(
                    f"{gid}: --proteins {proteins} unusable ({how}), "
                    f"annotating without it\n")
                proteins, src = "", f"{src}-unusable"
            else:
                if prepared != proteins:
                    sys.stderr.write(
                        f"{gid}: reformatted proteins for bakta ({how}): "
                        f"{prepared}\n")
                    src = f"{src}+reformatted"
                proteins = prepared

        rows.append({
            "genome_id": gid,
            "input_fasta": path,
            "assembly_type": kind,
            "domain": domain or "unknown",
            "tool": tool,
            "genus": genus,
            "species_epithet": epithet,
            "strain": strain,
            "gram": m.get("gram", "") if kind == "isolate_consensus" else "",
            "plasmid": m.get("plasmid", "") if kind == "isolate_consensus" else "",
            "gtdb_species": normalize_taxon(gtdb_species),
            "proteins": proteins,
            "proteins_source": src,
            "out_dir": str(Path(args.annotate_dir) / gid),
        })

    fields = ["genome_id", "input_fasta", "assembly_type", "domain", "tool",
              "genus", "species_epithet", "strain", "gram", "plasmid",
              "gtdb_species", "proteins", "proteins_source", "out_dir"]
    # Empty fields are written as a placeholder, NOT as an empty string.
    # bash's `IFS=$'\t' read` collapses runs of tabs (tab is IFS-whitespace),
    # so a genuinely empty column in the middle of a row would shift every
    # later column left and silently corrupt the annotator's arguments.
    # The consuming loop converts NONE_TOKEN back to "".
    # lineterminator="\n" is essential: csv defaults to "\r\n", and with
    # newline="" that CRLF reaches the file verbatim. bash's `read` strips
    # the \n but NOT the \r, so the last column (out_dir) would carry a
    # trailing carriage return and annotation would create directories
    # literally named "<genome_id>\r" that nothing downstream can find.
    NONE_TOKEN = "__NONE__"
    with open(args.out_plan, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (str(r.get(k, "")).strip() or NONE_TOKEN)
                        for k in fields})

    uniq_missing = sorted(set(m for m in missing if m))
    with open(args.out_missing, "w") as fh:
        for m in uniq_missing:
            fh.write(m + "\n")

    n_bak = sum(1 for r in rows if r["tool"] == "bakta")
    n_pro = sum(1 for r in rows if r["tool"] == "prokka")
    n_prot = sum(1 for r in rows if r["proteins"])
    print(f"planned {len(rows)} genomes: {n_bak} bakta (Bacteria), "
          f"{n_pro} prokka (Archaea); {n_prot} with --proteins, "
          f"{len(rows)-n_prot} without")
    if uniq_missing:
        print(f"unresolved species ({len(uniq_missing)}): "
              + ", ".join(uniq_missing[:8])
              + (" ..." if len(uniq_missing) > 8 else ""))


if __name__ == "__main__":
    main()
'''


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

    # 14. Replicon breakdown — main assembly vs Plassembler-validated plasmids vs bins
    if not isolates.empty:
        try:
            x = np.arange(len(isolates))
            plas_count = pd.to_numeric(isolates.get("plasmid_count", 0), errors="coerce").fillna(0)
            total_len  = pd.to_numeric(isolates.get("genome_length_bp", 0), errors="coerce").fillna(0)
            chrom_len  = pd.to_numeric(isolates.get("chromosome_length_bp", 0), errors="coerce").fillna(0)
            # When Plassembler validated >=1 plasmid, the chromosome bar is the
            # largest contig and plasmid bar is the sum of other contigs.
            # When Plassembler validated 0 plasmids, show the FULL consensus
            # as the main bar (no plasmid layer) — matches the report cards.
            has_plas = plas_count > 0
            chrom = np.where(has_plas, chrom_len, total_len) / 1e6
            plas  = pd.to_numeric(isolates.get("plasmid_total_bp", 0), errors="coerce").fillna(0) / 1e6
            plas  = np.where(has_plas, plas, 0)
            # Sum bin lengths per sample
            if not bins.empty:
                sample_bin_len = bins.groupby("sample")["genome_length_bp"].apply(
                    lambda s: pd.to_numeric(s, errors="coerce").sum()) / 1e6
                bin_vals = np.array([sample_bin_len.get(s, 0) for s in isolates["sample"]])
            else:
                bin_vals = np.zeros(len(isolates))
            fig, ax = plt.subplots(figsize=(11, max(5, 0.5*len(isolates)+1)))
            ax.barh(x, chrom, color="#2c7fb8", label="Main assembly")
            ax.barh(x, plas,  left=chrom, color="#7b68ee",
                    label="Plassembler-validated plasmid(s)")
            ax.barh(x, bin_vals, left=chrom + plas, color=COL_BIN,
                    label="MetaBAT2 bins")
            ax.set_yticks(x)
            ax.set_yticklabels(isolates["sample"])
            ax.invert_yaxis()
            ax.set_xlabel("Assembled length (Mb)")
            ax.set_title("Per-sample replicon breakdown")
            ax.legend(loc="lower right")
            # Annotate counts on the right of each bar
            for i, (cl, pl, bn) in enumerate(zip(chrom, plas, bin_vals)):
                n_plas_i = int(isolates["plasmid_count"].iloc[i]) if "plasmid_count" in isolates.columns else 0
                n_bin_i  = int(bins[bins["sample"] == isolates["sample"].iloc[i]].shape[0]) if not bins.empty else 0
                main_lbl = ("chr" if n_plas_i > 0 else "consensus")
                ax.text(cl + pl + bn + 0.05, i,
                        f"{main_lbl} {cl:.2f} Mb · {n_plas_i} plasmid(s) · {n_bin_i} bin(s)",
                        fontsize=8, va="center", color="#444")
            plots.append(("Replicon breakdown per sample",
                          "For each isolate: the main Autocycler assembly (whole "
                          "consensus if no plasmid was validated; chromosome only "
                          "if Plassembler validated plasmid contigs in the "
                          "clustering), the sum of any Plassembler-validated "
                          "plasmid contigs, and the total length of MetaBAT2 bins "
                          "recovered from unmapped reads.",
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
            # Render plasmid block only when Plassembler validated a plasmid.
            # Without Plassembler support, multi-contig consensuses are shown
            # as a single main assembly (no plasmid section at all).
            if plas_list:
                plas_rows = "".join(
                    f"<li><code>{escape(n)}</code> &mdash; {fmt_bp(l)}</li>"
                    for n, l in plas_list
                )
                plas_html = f"<ul class='replicon-list'>{plas_rows}</ul>"
                show_plasmid_section = True
            else:
                plas_html = ""
                show_plasmid_section = False

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

            # Main assembly description: show different info depending on
            # whether Plassembler validated any plasmid clusters.
            total_len = row.get("genome_length_bp", 0)
            n_total_contigs = int(float(row.get("n_contigs", 0) or 0))
            if show_plasmid_section:
                main_desc = (
                    f"<p>Chromosome contig: {fmt_bp(chrom_bp)}"
                    f"{' ' + badge if badge else ''}</p>"
                )
            else:
                contig_label = "1 contig" if n_total_contigs == 1 else f"{n_total_contigs} contigs"
                main_desc = (
                    f"<p>Consensus: {fmt_bp(total_len)} across {contig_label}"
                    f"{' ' + badge if badge else ''}</p>"
                )

            plasmid_block_html = ""
            if show_plasmid_section:
                plasmid_block_html = (
                    f"<div class='replicon-section'>"
                    f"<div class='replicon-label plasmids'>Plasmid contigs ({n_plas})"
                    f" &middot; Plassembler-validated</div>"
                    f"{plas_html}"
                    f"</div>"
                )

            sample_cards.append(
                f"<div class='sample-card'>"
                f"<h3>{escape(s)}</h3>"
                f"<p class='sample-meta'>Expected: <em>{escape(exp_str)}</em></p>"
                f"<div class='replicon-section'>"
                f"<div class='replicon-label main-asm'>Main assembly (Autocycler consensus)</div>"
                f"{main_desc}"
                f"<p>Classification: {fmt_classification(classification)}</p>"
                f"</div>"
                f"{plasmid_block_html}"
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


# ---- s21 MIJAMP preprocess ----------------------------------------------

def stage_mijamp(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    """Run MIJAMP's preprocess script per sample using the methylation-tagged
    BAM from s02b_demux_bam and the corresponding Bakta-annotated .fna from
    s16_annotate. MIJAMP is expected to be git-cloned at cfg['mijamp_dir'].
    """
    name = "s21_mijamp"
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg["sge_h_vmem_cpu"],
        cfg["sge_pe_parallel"], cfg["email"], hold_jid=hold,
    )
    body += conda_block(cfg["conda_env_mijamp"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    threads = int(cfg.get("mijamp_threads", 20))
    mij_dir = cfg.get("mijamp_dir", "") or ""
    body += textwrap.dedent(
        f"""
        MIJAMP_DIR="{mij_dir}"
        if [ -z "$MIJAMP_DIR" ] || [ ! -f "$MIJAMP_DIR/scripts/preprocess" ]; then
            log_event WARN "MIJAMP installation not found at '$MIJAMP_DIR' (set mijamp_dir in pipeline_config.yaml). Skipping."
            exit 0
        fi
        chmod +x "$MIJAMP_DIR/scripts/preprocess" 2>/dev/null || true

        mkdir -p {layout.mijamp}

        """
    )
    for s in samples:
        out_dir = f"{layout.mijamp}/{s.sample}"
        fna     = f"{layout.annotate}/{s.sample}/{s.sample}.fna"
        pats = barcode_filename_patterns(s.barcode)
        bam_name_clauses = " -o ".join(f'-name "*{p}*.bam"' for p in pats)
        body += textwrap.dedent(
            f"""
            # ---- {s.sample} (barcode {s.barcode}, matches: {", ".join(pats)}) ----
            if [ -d {out_dir} ] && [ -n "$(ls -A {out_dir} 2>/dev/null)" ]; then
                log_event SKIP "{s.sample} (MIJAMP output dir not empty)"
            else
                # Locate the per-barcode BAM. dorado demux with --no-classify
                # writes one BAM per BC tag, with the barcode label in the name.
                # Accept both NB13 and barcode13 naming conventions. Don't
                # cap recursion depth — older runs may have MinKNOW-style
                # nested layout from before s02b's flatten step existed.
                bam_file=$(find {layout.bams} -type f \\( {bam_name_clauses} \\) 2>/dev/null | head -1)

                # MIJAMP needs the GENOME SEQUENCE, not the annotation, so the
                # annotator's .fna is a convenience rather than a requirement.
                # Fall back to the Autocycler consensus (identical sequence)
                # when annotation is absent, empty or unreadable — a crashed
                # Bakta run should not also cost us the methylation analysis.
                # Look for the genome, retrying briefly. On NFS a file written
                # moments earlier by another node can stay invisible while the
                # attribute cache is warm, so a single stat is not conclusive.
                genome=""
                gsrc=""
                for _try in 1 2 3; do
                    if [ -s {fna} ]; then
                        genome={fna}; gsrc="annotation"; break
                    fi
                    if [ -s {layout.consensus}/{s.sample}.fasta ]; then
                        genome={layout.consensus}/{s.sample}.fasta; gsrc="consensus"
                        if [ -e {fna} ]; then
                            log_event WARN "{s.sample}: {fna} exists but is EMPTY (partial/failed annotation) — using the consensus instead"
                        fi
                        break
                    fi
                    [ "$_try" -lt 3 ] && {{ ls {layout.annotate}/{s.sample}/ >/dev/null 2>&1; sleep 5; }}
                done

                if [ -z "$bam_file" ]; then
                    log_event WARN "{s.sample}: no demux BAM matching any of: {', '.join(pats)}, skipping"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                elif [ -z "$genome" ]; then
                    log_event WARN "{s.sample}: no usable genome fasta after 3 attempts — neither {fna} nor {layout.consensus}/{s.sample}.fasta is present and non-empty; skipping"
                    log_event WARN "{s.sample}: dir listing was: $(ls -l {layout.annotate}/{s.sample}/ 2>&1 | tr '\\n' ' ' | cut -c1-300)"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    log_event OK "{s.sample}: MIJAMP preprocess (bam=$bam_file, genome=$genome [$gsrc])"
                    mkdir -p {out_dir}
                    set +e
                    "$MIJAMP_DIR/scripts/preprocess" \\
                        -b "$bam_file" \\
                        -g "$genome" \\
                        -t {threads} \\
                        -o {out_dir}
                    rc=$?
                    set -e
                    if [ $rc -ne 0 ]; then
                        log_event WARN "{s.sample}: MIJAMP rc=$rc"
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    else
                        log_event OK "{s.sample}: MIJAMP done"
                    fi
                fi
            fi

            """
        )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
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
        layout.annotate, layout.bams, layout.mapping, layout.meta_flye,
        layout.bins, layout.checkm2, layout.gtdbtk, layout.reports,
        layout.mijamp,
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
    # Parallel branch: dorado demux to BAM (preserves MM/ML methylation tags).
    # Runs in parallel with s02_demux off the same s01_basecall output, then
    # feeds s21_mijamp at the end.
    n02b = stage_demux_bam(cfg, layout, hold=n01)
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
    # GTDB-Tk and CheckM2 both take the raw assemblies (consensuses + bins)
    # and neither needs the other, so they fan out in parallel off s15_bin.
    n16 = stage_classify(cfg, layout, samples, hold=n15)
    n17 = stage_checkm2(cfg, layout, hold=n15)
    # Annotation now genuinely depends on classification: GTDB-Tk's assignment
    # picks BOTH the annotator (Bacteria->Bakta, Archaea->Prokka) and the
    # RefSeq type-strain proteome passed to --proteins.
    n18 = stage_annotate(cfg, layout, samples, hold=n16)
    # Aggregate reads annotation features + CheckM2 quality + GTDB-Tk taxonomy.
    n19 = stage_aggregate(cfg, layout, samples, hold=f"{n17},{n18}")
    n20 = stage_report(cfg, layout, hold=n19)
    # MIJAMP needs both the demux BAMs (s02b) and the Bakta .fna (s16, which
    # is upstream of s20). Hold on both so SGE waits for the slower of the two.
    n21 = stage_mijamp(cfg, layout, samples, hold=f"{n20},{n02b}")

    ordered = [
        [n01], [n02], [n02b], [n03], [n04], [n05],
        n06_list,
        [n07], [n08], [n09], [n10], [n11], [n12],
        [n13], [n14], [n15], [n16], [n17], [n18], [n19], [n20], [n21],
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
