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
  s01_basecall       dorado basecaller (POD5 -> basecalls.bam); with
                     basecall_split_jobs > 1 an SGE array, one GPU job per chunk
                     of pod5 files (each on its own GPU node)
  s01b_merge_bam     (basecall_split_jobs > 1) check the chunk BAMs and
                     concatenate them into basecalls.bam before demux
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
  s12b_polish        (polish_enabled) dorado polish --bacteria on each consensus,
                     before mapping and annotation (polish_tools.py; see README)
  s13_annotate       Bakta per sample using samplesheet metadata
  s14_classify       GTDB-Tk classify_wf across all Bakta .fna outputs
  s02b_demux_bam     parallel demux that emits sorted BAMs (retains MM/ML
                     methylation tags from basecalling)
  s18b_sv            structural heterogeneity within each sample: candidate
                     SVs (known genome versions and/or Sniffles2) re-genotyped
                     against explicit alleles with spanning reads or, for
                     events longer than the reads, at their breakpoint
                     junctions (sv_genotype.py; see README)
  s18c_snv           (snv_enabled) SNPs and indels between samples and against
                     reference genomes: reads of each sample mapped to its own
                     consensus, to near-identical samples and to snv_references,
                     validated caller configuration (snv_tools.py; see README)
  s18d_minor         (minor_enabled) minority SNVs within each sample against the
                     sample's own error model, with linkage into haplotype groups
                     (minor_variants.py; see README)
  s18e_phava         (phava_enabled) phase variation: invertible DNA segments
                     between inverted repeats, found and quantified with PhaVa on
                     each consensus and its reads, with gene context
                     (phava_tools.py; see README)
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
import hashlib
import os
import shutil
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
    array_tasks: int | None = None,
) -> str:
    # An array job writes one log per task; SGE substitutes $TASK_ID in -o.
    log_name = f"{job_name}.$TASK_ID" if array_tasks else job_name
    lines = [
        "#!/bin/bash",
        f"#$ -N {job_name}",
        f"#$ -o {stdout_dir}/{log_name}.log",
        "#$ -j y",
        "#$ -cwd",
        f"#$ -l h_rt={h_rt}",
        f"#$ -l h_vmem={h_vmem}",
        f"#$ -pe parallel {pe_parallel}",
        "#$ -m a",
        f"#$ -M {email}",
    ]
    if array_tasks:
        lines.append(f"#$ -t 1-{array_tasks}")
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
        # SGE appends every run of a stage to the same -o file, so note where
        # this run starts; log_warn_errors_from scans only the lines after it.
        STAGE_LOG_START_LINE=$( {{ wc -l < "$STAGE_LOG"; }} 2>/dev/null || echo 0 )
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
            # For the stage's own SGE log, skip the lines left by earlier runs
            # (STAGE_LOG_START_LINE, recorded at stage start); other files are
            # scanned in full. Reported line numbers stay absolute.
            local skip=0
            if [ "$f" = "$STAGE_LOG" ]; then skip="${{STAGE_LOG_START_LINE:-0}}"; fi
            skip="${{skip//[!0-9]/}}"; skip="${{skip:-0}}"
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
            tail -n +"$((skip + 1))" "$f" 2>/dev/null \\
                | grep -nE '(ERROR:|Error:|FATAL:|Fatal:|FAIL:|panic:|Segmentation fault|core dumped|^Killed$)' \\
                | awk -F: -v s="$skip" 'BEGIN {{ OFS = ":" }} {{ $1 = $1 + s; print }}' \\
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
    sv: Path
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
            sv=root / "sv",
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


# ---- s01 basecall split across GPU nodes + s01b merge ---------------------
#
# dorado cannot spread one job across nodes, so with basecall_split_jobs > 1 the
# pod5 files are divided into that many chunks at generate time and s01_basecall
# becomes an SGE array: task i basecalls chunk i on its own GPU node (the gpu
# complex is 1 per host here) into basecalled/chunks/chunk_<i>.bam. s01b_merge_bam
# then checks every chunk and concatenates them into basecalled/basecalls.bam, so
# s02_demux, s02b_demux_bam and everything after them are unchanged.

# Slowest card class measured on this cluster for sup v5.2 + 4mC_5mC + 6mA: 2.89-2.93
# GB of pod5 per hour on 4x TITAN Xp nodes (basecall_split_validation/), rounded down.
# Used only to warn when a chunk might not finish inside basecall_split_h_rt.
BASECALL_SLOWEST_GB_PER_HOUR = 2.8


def basecall_split_jobs(cfg: dict) -> int:
    """Number of GPU jobs s01_basecall is split into (1 = one job, the default)."""
    v = cfg.get("basecall_split_jobs", 1)
    if isinstance(v, bool) or not str(v).strip().isdigit() or int(v) < 1:
        raise SystemExit(f"basecall_split_jobs must be a whole number >= 1, got {v!r}")
    return int(v)


def _balance(files: List[tuple], n_jobs: int) -> List[List[tuple]]:
    """Largest file first, each onto the chunk with the fewest bytes so far (ties:
    lowest chunk number), so chunk sizes differ by at most one file."""
    n = min(n_jobs, len(files))
    chunks: List[List[tuple]] = [[] for _ in range(n)]
    load = [0] * n
    for name, size in sorted(files, key=lambda f: (-f[1], f[0])):
        i = min(range(n), key=lambda k: (load[k], k))
        chunks[i].append((name, size))
        load[i] += size
    for ch in chunks:
        ch.sort()
    return chunks


def _chunk_width(n: int) -> int:
    return max(2, len(str(n)))


def plan_basecall_chunks(pod5_dir: Path, n_jobs: int, previous: Path) -> tuple:
    """Assign the pod5 files dorado would read to chunks; returns (chunks, n_kept).

    dorado reads the *.pod5 files directly inside the folder (no --recursive), so
    those are the files assigned, to at most n_jobs balanced chunks. The plan is
    written to submit_scripts/s01_basecall_chunks.tsv, and s01b_merge_bam refuses
    to merge if the folder no longer matches it.

    Re-generating into the same output folder (for example after the sequencing run
    wrote more pod5 files) keeps the previous plan's chunks unchanged when all their
    files are still there at the same size, so finished chunks are not basecalled
    again; only the files the plan lacks are assigned, to up to n_jobs new chunks
    after the old ones. Delete the plan file to plan from scratch.
    """
    files = {p.name: p.stat().st_size for p in pod5_dir.iterdir()
             if p.name.endswith(".pod5") and p.is_file()}
    if not files:
        raise SystemExit("basecall_split_jobs > 1 plans the chunks from the pod5 files, "
                         f"but there are none directly inside {pod5_dir}")
    old: dict = {}
    if previous.is_file():
        with previous.open() as fh:
            next(fh, None)
            for line in fh:
                ch, name, size = line.rstrip("\n").split("\t")
                old.setdefault(ch, []).append((name, int(size)))
    if old:
        ids = sorted(old)
        planned = {name for ch in old.values() for name, _ in ch}
        new = sorted((name, size) for name, size in files.items() if name not in planned)
        extra = _balance(new, n_jobs) if new else []
        total = len(ids) + len(extra)
        if (ids == [f"{i:0{len(ids[0])}d}" for i in range(1, len(ids) + 1)]
                and len(ids[0]) == _chunk_width(total)
                and all(files.get(name) == size for ch in old.values() for name, size in ch)):
            return [sorted(old[i]) for i in ids] + extra, len(ids)
        print(f"{previous.name}: the previous chunk plan no longer matches {pod5_dir}; "
              "planning afresh, so finished chunks will be basecalled again")
    return _balance(sorted(files.items()), n_jobs), 0


def basecall_settings_key(cfg: dict) -> str:
    """Digest of the settings that change what dorado writes for a chunk; a chunk
    finished under other settings is basecalled again."""
    parts = (cfg["dorado_binary"], cfg["basecall_model"], ",".join(cfg["mod_models"]), cfg["kit_name"])
    return hashlib.sha256("\n".join(str(p) for p in parts).encode()).hexdigest()[:16]


def _h_rt_hours(h_rt: str):
    try:
        if ":" in h_rt:
            h, m, sec = (h_rt.split(":") + ["0", "0"])[:3]
            return int(h) + int(m) / 60 + int(sec) / 3600
        return int(h_rt) / 3600
    except ValueError:
        return None


def report_basecall_plan(cfg: dict, chunks: List[List[tuple]], requested: int, n_kept: int) -> None:
    sizes = [sum(size for _, size in ch) for ch in chunks]
    n_files = sum(len(ch) for ch in chunks)
    print(f"s01_basecall: {len(chunks)} GPU jobs over {n_files} pod5 files "
          f"({sum(sizes) / 1e9:.1f} GB); chunks of {min(sizes) / 1e9:.1f}-{max(sizes) / 1e9:.1f} GB")
    if n_kept == len(chunks):
        print("  same chunk plan as before (no pod5 files added)")
    elif n_kept:
        print(f"  kept chunks 1-{n_kept} of the previous plan; {len(chunks) - n_kept} new chunk(s) "
              "for pod5 files added since")
    elif len(chunks) < requested:
        print(f"  basecall_split_jobs is {requested}, but there are only {n_files} pod5 files")
    h_rt = str(cfg.get("basecall_split_h_rt") or cfg["sge_h_rt_gpu"])
    hours = _h_rt_hours(h_rt)
    worst = max(sizes) / 1e9 / BASECALL_SLOWEST_GB_PER_HOUR
    if hours is not None and worst > hours:
        print(f"WARNING: the largest chunk ({max(sizes) / 1e9:.1f} GB) can take {worst:.1f} h on the "
              f"slowest GPU class measured here ({BASECALL_SLOWEST_GB_PER_HOUR:g} GB/h), more than "
              f"basecall_split_h_rt ({h_rt}). Raise basecall_split_jobs or basecall_split_h_rt.",
              file=sys.stderr)


def stage_basecall_split(cfg: dict, layout: Layout, chunks: List[List[tuple]]) -> str:
    name = "s01_basecall"
    n = len(chunks)
    width = _chunk_width(n)
    manifest = layout.submit_dir / "s01_basecall_chunks.tsv"
    with manifest.open("w") as fh:
        fh.write("chunk\tpod5_file\tbytes\n")
        for i, ch in enumerate(chunks, 1):
            for fname, size in ch:
                fh.write(f"{i:0{width}d}\t{fname}\t{size}\n")
    mod_models = ",".join(cfg["mod_models"])
    chunk_dir = layout.basecalled / "chunks"
    body = sge_header(
        name, layout.stdout_dir,
        str(cfg.get("basecall_split_h_rt") or cfg["sge_h_rt_gpu"]),
        str(cfg.get("basecall_split_h_vmem") or cfg["sge_h_vmem_gpu"]),
        int(cfg.get("basecall_split_slots") or cfg["sge_pe_parallel"]),
        cfg["email"], gpu=True,
        gpumem=str(cfg.get("basecall_split_gpumem") or cfg["sge_gpumem"]),
        array_tasks=n,
    )
    body += conda_block(cfg["conda_env_basecall"]) + "\n"
    body += cuda_block(cfg) + "\n\n"
    # Each task logs under its own stage name (s01_basecall.<task>), so
    # pipeline_status.sh lists the chunks that are still running.
    body += log_setup(name + ".${SGE_TASK_ID}", layout) + "\n\n"
    body += textwrap.dedent(
        f"""
        TASK="${{SGE_TASK_ID:-}}"
        case "$TASK" in
            ''|*[!0-9]*)
                log_event WARN "s01_basecall.sh is an SGE array job; submit it with qsub (SGE_TASK_ID='$TASK')"
                exit 1 ;;
        esac
        CH=$(printf '%0{width}d' "$TASK")
        MANIFEST={manifest}
        CHUNK_DIR={chunk_dir}
        CHUNK_BAM="$CHUNK_DIR/chunk_$CH.bam"
        CHUNK_DONE="$CHUNK_DIR/chunk_$CH.done"
        LINKS="$CHUNK_DIR/pod5_$CH"
        # A finished chunk counts only for the same files and basecalling settings.
        SETTINGS_KEY={basecall_settings_key(cfg)}
        CH_LINES=$(awk -F'\\t' -v c="$CH" 'NR > 1 && $1 == c' "$MANIFEST")
        CHUNK_KEY=$(printf '%s\\n%s\\n' "$CH_LINES" "$SETTINGS_KEY" | sha256sum | cut -c1-64)
        mkdir -p "$CHUNK_DIR"

        if [ -z "$CH_LINES" ]; then
            log_event WARN "chunk $CH has no pod5 files in $MANIFEST"
            exit 1
        elif [ -s {layout.bam_path} ]; then
            log_event SKIP "basecalls.bam already exists ($(stat -c %s {layout.bam_path}) bytes); chunk $CH not needed"
        elif [ -s "$CHUNK_BAM" ] && [ -s "$CHUNK_DONE" ] && grep -qx "chunk_key=$CHUNK_KEY" "$CHUNK_DONE"; then
            log_event SKIP "chunk $CH already basecalled ($(stat -c %s "$CHUNK_BAM") bytes)"
        else
            rm -f "$CHUNK_DONE" "$CHUNK_BAM" "$CHUNK_BAM.partial"
            rm -rf "$LINKS"
            mkdir -p "$LINKS"
            n_files=0; n_bytes=0; n_bad=0
            while IFS=$'\\t' read -r c f sz; do
                src="{layout.pod5_dir}/$f"
                if [ ! -f "$src" ]; then
                    log_event WARN "chunk $CH: $f is no longer in {layout.pod5_dir}"
                    n_bad=$((n_bad+1)); continue
                fi
                now=$(stat -L -c %s "$src")
                if [ "$now" != "$sz" ]; then
                    log_event WARN "chunk $CH: $f is $now bytes but was $sz when the pipeline was generated"
                    n_bad=$((n_bad+1)); continue
                fi
                ln -s "$src" "$LINKS/$f"
                n_files=$((n_files+1)); n_bytes=$((n_bytes+sz))
            done <<< "$CH_LINES"
            if [ "$n_bad" -gt 0 ]; then
                log_event WARN "chunk $CH not basecalled: $n_bad pod5 file(s) missing or changed; re-run generate_pipeline.py, then submit_all.sh"
                rm -rf "$LINKS"
                exit 1
            fi
            gpus=""
            if command -v nvidia-smi >/dev/null 2>&1; then
                gpus=$(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | paste -sd ';' -)
            fi
            log_event OK "chunk $CH: dorado basecaller on $n_files pod5 files ($n_bytes bytes) on $(hostname), GPUs: ${{gpus:-not reported}}"
            t0=$(date +%s)
            {cfg['dorado_binary']} basecaller \\
                {cfg['basecall_model']} \\
                "$LINKS" \\
                --modified-bases-models {mod_models} \\
                --kit-name {cfg['kit_name']} \\
                > "$CHUNK_BAM.partial"
            rc=$?
            t1=$(date +%s)
            if [ $rc -ne 0 ] || [ ! -s "$CHUNK_BAM.partial" ]; then
                log_event WARN "dorado basecaller exited with code $rc on chunk $CH"
                rm -f "$CHUNK_BAM.partial"
                rm -rf "$LINKS"
                exit $(( rc == 0 ? 1 : rc ))
            fi
            mv -f "$CHUNK_BAM.partial" "$CHUNK_BAM"
            {{
                echo "chunk=$CH"
                echo "chunk_key=$CHUNK_KEY"
                echo "pod5_files=$n_files"
                echo "pod5_bytes=$n_bytes"
                echo "bam_bytes=$(stat -c %s "$CHUNK_BAM")"
                echo "host=$(hostname)"
                echo "gpus=$gpus"
                echo "start_epoch=$t0"
                echo "seconds=$((t1 - t0))"
            }} > "$CHUNK_DONE.tmp" && mv -f "$CHUNK_DONE.tmp" "$CHUNK_DONE"
            rm -rf "$LINKS"
            log_event OK "chunk $CH basecalled in $((t1 - t0)) s ($(stat -c %s "$CHUNK_BAM") bytes)"
        fi

        log_warn_errors_from "$STAGE_LOG"
        """
    ).strip() + "\n"
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


def stage_merge_basecalls(cfg: dict, layout: Layout, chunks: List[List[tuple]], hold: str) -> str:
    name = "s01b_merge_bam"
    n = len(chunks)
    width = _chunk_width(n)
    threads = 4
    manifest = layout.submit_dir / "s01_basecall_chunks.tsv"
    chunk_dir = layout.basecalled / "chunks"
    keep = 1 if cfg.get("basecall_split_keep_chunks", False) else 0
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"],
        cfg.get("sge_h_vmem_merge_bam", "4G"), threads, cfg["email"], hold_jid=hold,
    )
    # Activating the nanopore env runs `R CMD javareconf`, and R's threaded OpenBLAS reserves
    # address space for every core of the node: under this stage's small h_vmem on a 128-core
    # node the activation then spins instead of finishing (>300 s, against 2 s with one thread).
    body += "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\n"
    body += conda_block(cfg["conda_env_demux"]) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    body += textwrap.dedent(
        f"""
        if ! check_samtools "{cfg['conda_env_demux']}"; then
            log_event WARN "aborting s01b - samtools >= 1.10 is required to merge the chunk BAMs"
            exit 1
        fi
        # quickcheck -u: dorado's BAMs are unaligned, and without -u quickcheck rejects
        # any BAM whose header lists no reference sequences.
        MANIFEST={manifest}
        CHUNK_DIR={chunk_dir}
        OUT={layout.bam_path}
        N={n}
        SETTINGS_KEY={basecall_settings_key(cfg)}
        KEEP_CHUNKS={keep}
        mkdir -p "$CHUNK_DIR"

        if [ -s "$OUT" ]; then
            log_event SKIP "basecalls.bam already exists ($(stat -c %s "$OUT") bytes)"
        else
            n_bad=0
            # 1. The pod5 folder still holds exactly the files the chunks were planned from.
            find -L {layout.pod5_dir} -maxdepth 1 -name '*.pod5' -type f -printf '%f\\n' | LC_ALL=C sort > "$CHUNK_DIR/pod5_now.txt"
            awk -F'\\t' 'NR > 1 {{ print $2 }}' "$MANIFEST" | LC_ALL=C sort > "$CHUNK_DIR/pod5_planned.txt"
            n_new=$(LC_ALL=C comm -23 "$CHUNK_DIR/pod5_now.txt" "$CHUNK_DIR/pod5_planned.txt" | wc -l)
            n_gone=$(LC_ALL=C comm -13 "$CHUNK_DIR/pod5_now.txt" "$CHUNK_DIR/pod5_planned.txt" | wc -l)
            if [ "$n_new" -gt 0 ]; then
                first_new=$(LC_ALL=C comm -23 "$CHUNK_DIR/pod5_now.txt" "$CHUNK_DIR/pod5_planned.txt" | head -1)
                running=""
                ls {layout.pod5_dir}/../final_summary_*.txt >/dev/null 2>&1 || running=" (no final_summary_*.txt next to the pod5 folder: the sequencing run may still be writing files)"
                log_event WARN "$n_new pod5 file(s) in {layout.pod5_dir} are not in the chunk plan, e.g. $first_new$running; their reads would be missing. Re-run generate_pipeline.py with the same arguments (it keeps the finished chunks and adds chunks for the new files), then submit_all.sh"
                n_bad=$((n_bad+1))
            fi
            if [ "$n_gone" -gt 0 ]; then
                log_event WARN "$n_gone pod5 file(s) of the chunk plan are no longer in {layout.pod5_dir}"
                n_bad=$((n_bad+1))
            fi

            # 2. Every chunk finished with the current plan and settings, and its BAM is intact.
            printf 'chunk\\tpod5_files\\tpod5_bytes\\treads\\tbam_bytes\\thost\\tgpus\\tseconds\\tgb_per_hour\\n' > "$CHUNK_DIR/summary.tsv.tmp"
            done_field() {{ awk -v k="$1" 'index($0, k "=") == 1 {{ print substr($0, length(k) + 2); exit }}' "$2"; }}
            total=0
            inputs=()
            for i in $(seq 1 "$N"); do
                CH=$(printf '%0{width}d' "$i")
                bam="$CHUNK_DIR/chunk_$CH.bam"
                done_f="$CHUNK_DIR/chunk_$CH.done"
                ch_lines=$(awk -F'\\t' -v c="$CH" 'NR > 1 && $1 == c' "$MANIFEST")
                key=$(printf '%s\\n%s\\n' "$ch_lines" "$SETTINGS_KEY" | sha256sum | cut -c1-64)
                if [ ! -s "$bam" ] || [ ! -s "$done_f" ]; then
                    log_event WARN "chunk $CH is not finished (see {layout.stdout_dir}/s01_basecall.$i.log)"
                    n_bad=$((n_bad+1)); continue
                fi
                if ! grep -qx "chunk_key=$key" "$done_f"; then
                    log_event WARN "chunk $CH was basecalled for another chunk plan or other basecalling settings; re-run submit_all.sh to redo it"
                    rm -f "$done_f"
                    n_bad=$((n_bad+1)); continue
                fi
                if ! samtools quickcheck -u "$bam"; then
                    log_event WARN "chunk $CH BAM is truncated or corrupt; moved to chunk_$CH.bam.bad, re-run submit_all.sh to basecall it again"
                    mv -f "$bam" "$bam.bad"
                    rm -f "$done_f"
                    n_bad=$((n_bad+1)); continue
                fi
                reads=$(samtools view -c -@ {threads} "$bam")
                if ! [ "$reads" -ge 0 ] 2>/dev/null; then
                    log_event WARN "could not count the reads of chunk $CH"
                    n_bad=$((n_bad+1)); continue
                fi
                # Empty fields become NA: tab is whitespace to `read`, so an empty
                # field would shift the ones after it when the summary is read back.
                secs=$(done_field seconds "$done_f"); secs=${{secs:-NA}}
                pbytes=$(done_field pod5_bytes "$done_f"); pbytes=${{pbytes:-NA}}
                nfiles=$(done_field pod5_files "$done_f"); nfiles=${{nfiles:-NA}}
                host=$(done_field host "$done_f"); host=${{host:-NA}}
                gpus=$(done_field gpus "$done_f"); gpus=${{gpus:-NA}}
                gbh=$(awk -v b="$pbytes" -v s="$secs" 'BEGIN {{ if (s + 0 > 0 && b + 0 > 0) printf "%.2f", b / 1e9 / (s / 3600); else print "NA" }}')
                printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "$CH" "$nfiles" "$pbytes" "$reads" \\
                    "$(stat -c %s "$bam")" "$host" "$gpus" "$secs" "$gbh" >> "$CHUNK_DIR/summary.tsv.tmp"
                total=$((total + reads))
                inputs+=("$bam")
            done
            if [ "$n_bad" -gt 0 ]; then
                log_event WARN "basecalls.bam not written: $n_bad problem(s) above"
                exit 1
            fi

            # 3. Header: @HD and @PG of chunk 1, the read groups of every chunk (an id must
            #    carry the same line in all of them), and @CO lines recording the split.
            HDR="$CHUNK_DIR/merged_header.sam"
            samtools view -H "${{inputs[0]}}" | awk '/^@HD/' > "$HDR.tmp"
            for b in "${{inputs[@]}}"; do samtools view -H "$b"; done | awk -F'\\t' '
                /^@RG/ {{
                    id = ""
                    for (k = 2; k <= NF; k++) if ($k ~ /^ID:/) id = $k
                    if (id in seen) {{ if (seen[id] != $0) bad = 1; next }}
                    seen[id] = $0
                    print
                }}
                END {{ exit bad }}' >> "$HDR.tmp"
            rg_rc=$?
            if [ "$rg_rc" -ne 0 ]; then
                log_event WARN "chunks disagree on the header line of a read group; basecalls.bam not written"
                rm -f "$HDR.tmp"
                exit 1
            fi
            samtools view -H "${{inputs[0]}}" | awk '/^@PG/' >> "$HDR.tmp"
            printf '@CO\\ts01_basecall ran as %s GPU jobs over the pod5 files in %s (chunk plan %s); the basecaller @PG line is that of chunk 1, every chunk ran the same command on links to its own pod5 files\\n' \\
                "$N" "{layout.pod5_dir}" "$MANIFEST" >> "$HDR.tmp"
            tail -n +2 "$CHUNK_DIR/summary.tsv.tmp" | while IFS=$'\\t' read -r ch nf pb rd bb host gpus secs gbh; do
                printf '@CO\\ts01_basecall chunk %s: %s pod5 files, %s bytes, %s reads, basecalled in %s s on %s (%s)\\n' \\
                    "$ch" "$nf" "$pb" "$rd" "$secs" "$host" "$gpus"
            done >> "$HDR.tmp"
            mv -f "$HDR.tmp" "$HDR"

            # 4. Concatenate in chunk order and check the result before it is put in place.
            log_event OK "concatenating $N chunk BAMs ($total reads)"
            samtools cat -h "$HDR" -o "$OUT.partial" "${{inputs[@]}}"
            rc=$?
            if [ "$rc" -ne 0 ] || ! samtools quickcheck -u "$OUT.partial"; then
                log_event WARN "samtools cat failed (rc=$rc); basecalls.bam not written"
                rm -f "$OUT.partial"
                exit 1
            fi
            n_out=$(samtools view -c -@ {threads} "$OUT.partial")
            if [ "$n_out" != "$total" ]; then
                log_event WARN "the merged BAM has $n_out records but the chunks have $total; basecalls.bam not written"
                rm -f "$OUT.partial"
                exit 1
            fi
            mv -f "$OUT.partial" "$OUT"
            mv -f "$CHUNK_DIR/summary.tsv.tmp" "$CHUNK_DIR/summary.tsv"
            log_event OK "basecalls.bam written: $total reads from $N chunks ($(stat -c %s "$OUT") bytes); per-chunk times in $CHUNK_DIR/summary.tsv"
            if [ "$KEEP_CHUNKS" = 1 ]; then
                log_event OK "kept the chunk BAMs in $CHUNK_DIR (basecall_split_keep_chunks: true)"
            else
                # The chunk BAMs duplicate basecalls.bam. Their .done records stay, and
                # the tasks skip on a re-run because basecalls.bam exists.
                rm -f "${{inputs[@]}}"
                log_event OK "removed the $N chunk BAMs (basecall_split_keep_chunks: false)"
            fi
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


# ---- s12b consensus polishing --------------------------------------------

# Defaults chosen by the validation rules (README "Small variants: validation").
POLISH_DEFAULT_DEPTH = 100
SNV_DEFAULT_CALLER = "INT"
SNV_DEFAULT_DEPTH = 25

POLISH_ENGINE = Path(__file__).resolve().parent / "polish_tools.py"


def stage_polish(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str | None:
    """Polish each Autocycler consensus with `dorado polish --bacteria` (s12b).

    Written only with polish_enabled: true; otherwise there is no s12b script and
    s13 holds on s12_collect as before. Reads: a deterministic read-name-hash
    subsample of fastq_filtered/<sample>.fastq at polish_depth-fold coverage,
    mapped with minimap2 -x lr:hq keeping each read's RG tag (dorado picks its
    polishing model from the basecall model in the @RG header). The change table
    comes from dorado's VCF and must turn the draft into dorado's consensus
    exactly; only then is consensus/<sample>.fasta replaced (same contig names
    and FASTA headers, length= updated). A consensus that downstream stages have
    already used (s13 mapping or annotation present) is never replaced unless
    polish_replace_existing is true. With polish_install: false (review mode)
    the polished genome and change table are written but the consensus is left
    unchanged.
    Outputs: polish/<sample>/{draft.fasta, polished.fasta, changes.tsv,
    variants.vcf, info.json, logs}.
    """
    if not bool(cfg.get("polish_enabled", False)):
        return None
    name = "s12b_polish"
    threads = int(cfg.get("polish_threads", 16))
    env = cfg.get("conda_env_snv", "snv")
    depth = int(cfg.get("polish_depth", POLISH_DEFAULT_DEPTH))
    device = str(cfg.get("polish_device", "cpu"))
    replace = bool(cfg.get("polish_replace_existing", False))
    install = bool(cfg.get("polish_install", True))
    max_change = float(cfg.get("polish_max_length_change", 0.001))
    model_name = Path(str(cfg["basecall_model"])).name
    models_dir = str(cfg.get("polish_models_dir") or Path(str(cfg["basecall_model"])).parent)
    gpu = device.startswith("cuda")
    problems = []
    if not POLISH_ENGINE.is_file():
        problems.append(f"{POLISH_ENGINE} is missing (it ships next to generate_pipeline.py)")
    if not 10 <= depth <= 1000:
        problems.append("polish_depth must be between 10 and 1000")
    if device != "cpu" and not gpu:
        problems.append(f"polish_device must be cpu or cuda:... (got {device!r})")
    if not 0 < max_change < 0.1:
        problems.append("polish_max_length_change must be in (0, 0.1)")
    if threads < 1:
        problems.append("polish_threads must be >= 1")
    if problems:
        raise SystemExit("pipeline_config.yaml (s12b_polish): " + "; ".join(problems))
    engine = layout.submit_dir / "_polish_tools.py"
    shutil.copyfile(POLISH_ENGINE, engine)
    engine.chmod(engine.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_gpu"] if gpu else cfg["sge_h_rt_cpu"],
        cfg.get("sge_h_vmem_polish", "6G"), threads, cfg["email"],
        gpu=gpu, gpumem=cfg.get("sge_gpumem", "20G"), hold_jid=hold,
    )
    body += conda_block(env) + "\n\n"
    if gpu:
        body += cuda_block(cfg) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    pol = layout.root / "polish"
    body += textwrap.dedent(
        f"""
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
        set +e
        set -o pipefail
        POLISH_DIR={pol}
        PT={engine}
        DORADO={cfg['dorado_binary']}
        MODELS={models_dir}
        mkdir -p "$POLISH_DIR"
        if ! python3 -c "import mappy, pysam" 2>/dev/null; then
            log_event WARN "the '{env}' conda env lacks mappy/pysam (needs python>=3.9, mappy, pysam, numpy, scipy, minimap2, samtools, bcftools); consensus genomes left unpolished"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            exit 0
        fi
        for tool in minimap2 samtools; do
            if ! command -v "$tool" >/dev/null 2>&1; then
                log_event WARN "$tool not found in the '{env}' env; consensus genomes left unpolished"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                exit 0
            fi
        done
        if [ ! -x "$DORADO" ]; then
            log_event WARN "dorado not found at $DORADO (dorado_binary); consensus genomes left unpolished"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            exit 0
        fi
        DORADO_VERSION=$("$DORADO" --version 2>&1 | tail -1)
        ls -d "$MODELS"/*polish* >/dev/null 2>&1 \\
            || log_event WARN "no polishing model under $MODELS (polish_models_dir); dorado will try to download one"
        log_event OK "dorado $DORADO_VERSION; polishing at {depth}x on {device} with {threads} threads{'' if install else '; review mode (polish_install: false): consensus genomes are not replaced'}"

        # One sample: everything happens in polish/<sample>.tmp, which replaces
        # polish/<sample> only when every step has succeeded.
        polish_sample() {{
            local S="$1" CONS="$2" READS="$3" W="$POLISH_DIR/$1.tmp" mode
            rm -rf "$W" && mkdir -p "$W" || return 1
            cp "$CONS" "$W/draft.fasta" || return 1
            python3 "$PT" stats --fastq "$READS" --out "$W/stats.json" >> "$W/polish.log" 2>&1 || return 1
            python3 "$PT" subsample --fastq "$READS" --stats "$W/stats.json" --draft "$W/draft.fasta" \\
                --depths {depth} --seed polish1 --prefix "$W/sub" >> "$W/polish.log" 2>&1 || return 1
            mode=$(python3 "$PT" rgmode --stats "$W/stats.json") || return 1
            if [ "$mode" = tags ]; then
                python3 "$PT" rgheader --stats "$W/stats.json" --sample "$S" --basecall-model "{model_name}" \\
                    --out "$W/rg.txt" >> "$W/polish.log" 2>&1 || return 1
                minimap2 -ax lr:hq -y -t {threads} "$W/draft.fasta" "$W/sub.d{depth}.fastq" 2>> "$W/minimap2.log" \\
                    | samtools view -b -o "$W/aln.u.bam" - || return 1
                {{ samtools view -H "$W/aln.u.bam" | awk '!/^@RG/'; cat "$W/rg.txt"; }} > "$W/aln.hdr" || return 1
                samtools reheader "$W/aln.hdr" "$W/aln.u.bam" | samtools sort -@ 4 -m 1G -o "$W/aln.bam" - || return 1
                rm -f "$W/aln.u.bam" "$W/aln.hdr"
            elif [ "$mode" = none ]; then
                minimap2 -ax lr:hq -R "@RG\\tID:$S\\tPL:ONT\\tDS:basecall_model={model_name}\\tSM:$S" -t {threads} \\
                    "$W/draft.fasta" "$W/sub.d{depth}.fastq" 2>> "$W/minimap2.log" \\
                    | samtools sort -@ 4 -m 1G -o "$W/aln.bam" - || return 1
            else
                echo "some reads carry RG tags and some do not; cannot give dorado one basecall model per read" >> "$W/polish.log"
                return 1
            fi
            samtools index "$W/aln.bam" || return 1
            "$DORADO" polish "$W/aln.bam" "$W/draft.fasta" --bacteria --models-directory "$MODELS" -x {device} \\
                -t {threads} -o "$W/dorado" --vcf --any-bam > "$W/dorado_polish.log" 2>&1 || return 1
            python3 "$PT" changes --draft "$W/draft.fasta" --polished "$W/dorado/consensus.fasta" --bam "$W/aln.bam" \\
                --vcf "$W/dorado/variants.vcf" --out "$W/changes.tsv" >> "$W/polish.log" 2>&1 || return 1
            python3 "$PT" finalize --sample "$S" --draft "$W/draft.fasta" --polished "$W/dorado/consensus.fasta" \\
                --changes "$W/changes.tsv" --stats "$W/stats.json" --subsample "$W/sub.subsample.json" \\
                --depth {depth} --max-length-change {max_change} --out-polished "$W/polished.fasta" \\
                --out-info "$W/info.json" --dorado-version "$DORADO_VERSION" >> "$W/polish.log" 2>&1 || return 1
            cp "$W/dorado/variants.vcf" "$W/variants.vcf" || return 1
            rm -rf "$W/dorado" "$W/sub.d{depth}.fastq" "$W/aln.bam" "$W/aln.bam.bai"
            rm -rf "$POLISH_DIR/$S" && mv "$W" "$POLISH_DIR/$S"
        }}
        """
    )
    for s in samples:
        cons = f"{layout.consensus}/{s.sample}.fasta"
        reads = f"{layout.fastq_filtered}/{s.sample}.fastq"
        if not install:
            body += textwrap.dedent(
                f"""
                # ---- {s.sample} (review mode) ----
                if [ ! -s {cons} ]; then
                    log_event WARN "{s.sample}: consensus missing; not polished"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                elif python3 "$PT" is-installed --as-draft --info "$POLISH_DIR/{s.sample}/info.json" --consensus {cons}; then
                    log_event SKIP "{s.sample} (polishing results for this consensus already in $POLISH_DIR/{s.sample})"
                elif [ ! -s {reads} ]; then
                    log_event WARN "{s.sample}: filtered reads missing; not polished"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    t_s=$(date +%s)
                    if polish_sample "{s.sample}" {cons} {reads}; then
                        log_event OK "{s.sample} polished at {depth}x for review (consensus unchanged): $(python3 "$PT" describe --info "$POLISH_DIR/{s.sample}/info.json") in $(( $(date +%s) - t_s )) s"
                    else
                        log_event WARN "{s.sample}: polishing failed; logs in $POLISH_DIR/{s.sample}.tmp: $(tail -n 3 "$POLISH_DIR/{s.sample}.tmp/polish.log" "$POLISH_DIR/{s.sample}.tmp/dorado_polish.log" 2>/dev/null | tr '\\n' ' ' | cut -c1-400)"
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    fi
                fi
                """
            )
            continue
        used = (f'[ -e {layout.mapping}/{s.sample}/.mapping_done ] || [ -e {layout.mapping}/{s.sample}/mapped.bam ] '
                f'|| [ -n "$(ls -A {layout.annotate}/{s.sample} 2>/dev/null)" ]')
        body += textwrap.dedent(
            f"""
            # ---- {s.sample} ----
            if [ ! -s {cons} ]; then
                log_event WARN "{s.sample}: consensus missing; not polished"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            elif python3 "$PT" is-installed --info "$POLISH_DIR/{s.sample}/info.json" --consensus {cons}; then
                log_event SKIP "{s.sample} (consensus is already the polished genome; delete $POLISH_DIR/{s.sample} and restore the draft to re-run)"
            elif [ "{str(replace).lower()}" != true ] && {{ {used}; }}; then
                log_event WARN "{s.sample}: consensus already used downstream (s13 mapping or annotation exist), so it was not replaced. To polish it, move {layout.mapping}/{s.sample} and {layout.annotate}/{s.sample} aside (or set polish_replace_existing: true) and re-run from s12b"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            elif [ ! -s {reads} ]; then
                log_event WARN "{s.sample}: filtered reads missing; not polished"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                t_s=$(date +%s)
                if ! polish_sample "{s.sample}" {cons} {reads}; then
                    log_event WARN "{s.sample}: polishing failed, consensus left unpolished; logs in $POLISH_DIR/{s.sample}.tmp: $(tail -n 3 "$POLISH_DIR/{s.sample}.tmp/polish.log" "$POLISH_DIR/{s.sample}.tmp/dorado_polish.log" 2>/dev/null | tr '\\n' ' ' | cut -c1-400)"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                elif ! python3 "$PT" install --info "$POLISH_DIR/{s.sample}/info.json" \\
                        --polished "$POLISH_DIR/{s.sample}/polished.fasta" --consensus {cons} >> "$POLISH_DIR/{s.sample}/polish.log" 2>&1; then
                    log_event WARN "{s.sample}: polished genome could not be installed as {cons}: $(tail -n 2 "$POLISH_DIR/{s.sample}/polish.log" | tr '\\n' ' ' | cut -c1-300)"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                else
                    log_event OK "{s.sample} polished at {depth}x: $(python3 "$PT" describe --info "$POLISH_DIR/{s.sample}/info.json") in $(( $(date +%s) - t_s )) s; draft kept as $POLISH_DIR/{s.sample}/draft.fasta"
                fi
            fi
            """
        )
    # pipefail is on for the steps above; the final scan's grep finds nothing on a clean
    # run and must not turn that into the script's exit status.
    body += '\nset +o pipefail\nlog_warn_errors_from "$STAGE_LOG"\n'
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


# ---- s18b structural heterogeneity ---------------------------------------

SV_ENGINE = Path(__file__).resolve().parent / "sv_genotype.py"


def stage_sv(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    """Structural heterogeneity inside each sample (two or more versions of
    the same genome in one culture), with sv_genotype.py.

    Candidates: differences between the genome versions listed in
    sv_variant_genomes (e.g. a strain and its deletion mutant) and, with
    sv_discovery, Sniffles2 calls on each sample's s13 mapping. Every
    candidate is re-genotyped against explicit allele sequences: a read counts
    for an allele only if one alignment spans that allele's whole window, and
    all alleles share one window length, so fractions are not biased by allele
    length. Events longer than the reads are genotyped at their breakpoint
    junctions instead (sv_junction_mode). Each sample gets an empirical
    detection floor from size-matched fake events at unique positions. Every
    locus is typed in every sample whose consensus contains its flanks.
    Outputs: sv/events.tsv, sv/null_loci.tsv, sv/summary.json,
    sv/sample_summary.tsv, sv/<sample>/ (alleles, read lists, Sniffles2 VCF).
    """
    name = "s18b_sv"
    threads = int(cfg.get("sv_threads", 8))
    env = cfg.get("conda_env_sv", "sv")
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg.get("sge_h_vmem_sv", "4G"),
        threads, cfg["email"], hold_jid=hold,
    )
    body += conda_block(env) + "\n\n"
    body += log_setup(name, layout) + "\n\n"

    enabled = bool(cfg.get("sv_enabled", True))
    discovery = bool(cfg.get("sv_discovery", True))
    variants = [os.path.expandvars(os.path.expanduser(str(v))) for v in (cfg.get("sv_variant_genomes") or [])]
    jmode = str(cfg.get("sv_junction_mode", "auto"))
    min_len = int(cfg.get("sv_min_len", 500))
    anchor = int(cfg.get("sv_anchor", 1000))
    min_alt = int(cfg.get("sv_min_alt_reads", 5))
    min_frac = float(cfg.get("sv_min_fraction", 0.01))
    min_span = int(cfg.get("sv_min_spanning", 20))
    n_null = int(cfg.get("sv_null_loci", 20))
    max_reads = int(cfg.get("sv_max_reads", 3000))
    max_span = int(cfg.get("sv_max_span_window", 20000))
    engine = layout.submit_dir / "_sv_genotype.py"
    if enabled:
        problems = []
        if not SV_ENGINE.is_file():
            problems.append(f"{SV_ENGINE} is missing (it ships next to generate_pipeline.py)")
        missing = [v for v in variants if not Path(v).is_file()]
        if missing:
            problems.append("sv_variant_genomes not found: " + ", ".join(missing))
        if len(variants) == 1:
            problems.append("sv_variant_genomes needs two or more genome versions (or none)")
        if jmode not in ("auto", "always", "never"):
            problems.append(f"sv_junction_mode must be auto, always or never (got {jmode!r})")
        if (not 0 < min_frac < 0.5 or min_len < 50 or anchor < 200 or min_alt < 1 or min_span < 1 or n_null < 0
                or max_span < 2 * anchor):
            problems.append("sv_min_fraction must be in (0, 0.5); sv_min_len >= 50; sv_anchor >= 200; "
                            "sv_min_alt_reads, sv_min_spanning >= 1; sv_null_loci >= 0; sv_max_span_window >= 2 x sv_anchor")
        if problems:
            raise SystemExit("pipeline_config.yaml (s18b_sv): " + "; ".join(problems))
        shutil.copyfile(SV_ENGINE, engine)
        engine.chmod(engine.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    body += textwrap.dedent(
        f"""
        # OpenBLAS (pulled in by numpy) reserves memory per core at import; under
        # h_vmem on many-core nodes that alone can exceed the limit.
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
        SV_DIR={layout.sv}
        mkdir -p "$SV_DIR"
        """
    )
    if not enabled:
        body += 'log_event SKIP "structural-heterogeneity stage disabled (sv_enabled: false)"\nexit 0\n'
        write_script(layout.submit_dir / f"{name}.sh", body)
        return name
    if not variants and not discovery:
        body += ('log_event SKIP "no candidate source: sv_variant_genomes is empty and sv_discovery is false"\n'
                 'exit 0\n')
        write_script(layout.submit_dir / f"{name}.sh", body)
        return name

    vg_arg = ("--variant-genomes " + " ".join(f'"{v}"' for v in variants)) if variants else ""
    opts = (f"--min-len {min_len} --anchor {anchor} --min-alt-reads {min_alt} --min-fraction {min_frac} "
            f"--min-spanning {min_span} --null-loci {n_null} --max-reads {max_reads} --junction-mode {jmode} "
            f"--max-span-window {max_span}")
    body += textwrap.dedent(
        f"""
        if [ -s "$SV_DIR/events.tsv" ] && [ -f "$SV_DIR/.complete" ]; then
            log_event SKIP "results already in $SV_DIR (delete that directory to re-run)"
            exit 0
        fi
        if ! python3 -c "import mappy, pysam" 2>/dev/null; then
            log_event WARN "the '{env}' conda env lacks mappy/pysam (needs python>=3.9, mappy, pysam, sniffles>=2.2, samtools); skipping"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            exit 0
        fi
        DISC="{'--discovery' if discovery else ''}"
        if [ -n "$DISC" ] && ! command -v sniffles >/dev/null 2>&1; then
            log_event WARN "sniffles not found in the '{env}' env; genotyping the supplied genome versions only"
            DISC=""
        fi
        if [ -z "$DISC" ] && [ -z '{vg_arg}' ]; then
            log_event WARN "no candidate source left (no sniffles, no sv_variant_genomes); skipping"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            exit 0
        fi
        : > "$SV_DIR/samples.tsv.tmp"
        """
    )
    for s in samples:
        cons = f"{layout.consensus}/{s.sample}.fasta"
        bam = f"{layout.mapping}/{s.sample}/mapped.bam"
        adir = f"{layout.annotate}/{s.sample}"
        body += textwrap.dedent(
            f"""
            # ---- {s.sample} ----
            gff=""
            for g in {adir}/{s.sample}.gff3 {adir}/{s.sample}.gff; do
                if [ -s "$g" ]; then gff="$g"; break; fi
            done
            fna={adir}/{s.sample}.fna
            [ -s "$fna" ] || fna=""
            if [ ! -s {cons} ] || [ ! -s {bam} ] || [ ! -s {bam}.bai ]; then
                log_event WARN "{s.sample}: consensus or indexed s13 mapping (mapped.bam + .bai) missing; not genotyped"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "{s.sample}" {cons} {bam} "$gff" "$fna" >> "$SV_DIR/samples.tsv.tmp"
                [ -n "$gff" ] || log_event WARN "{s.sample}: no annotation GFF; genes in variable segments will not be listed"
            fi
            """
        )
    body += textwrap.dedent(
        f"""
        mv "$SV_DIR/samples.tsv.tmp" "$SV_DIR/samples.tsv"
        n_s=$(wc -l < "$SV_DIR/samples.tsv")
        if [ "$n_s" -eq 0 ]; then
            log_event WARN "no sample has a consensus + indexed mapping; nothing to genotype"
            exit 0
        fi
        log_event OK "genotyping structural heterogeneity in $n_s sample(s): {len(variants)} supplied genome version(s), discovery $([ -n "$DISC" ] && echo on || echo off), junction mode {jmode}"
        set +e
        python3 {engine} --samples "$SV_DIR/samples.tsv" --outdir "$SV_DIR" $DISC {vg_arg} \\
            --threads {threads} {opts} > "$SV_DIR/sv_genotype.log" 2>&1
        rc=$?
        set -e
        if [ $rc -ne 0 ] || [ ! -s "$SV_DIR/events.tsv" ]; then
            log_event WARN "sv_genotype rc=$rc; last lines: $(tail -3 "$SV_DIR/sv_genotype.log" | tr '\\n' ' ' | cut -c1-400)"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
        else
            tail -n +2 "$SV_DIR/sample_summary.tsv" | while IFS=$'\\t' read -r smp nloc nkn ndi mix dis sgl nsup ins rl und rep; do
                log_event OK "$smp: $nloc loci ($nkn known, $ndi discovered); mixture $mix, consensus_discrepancy $dis, single $sgl, not_supported $nsup, insufficient reads $ins, repeat-limited $rl"
                [ "$rep" = "-" ] || log_event OK "$smp: $(printf '%s' "$rep" | cut -c1-900)"
            done
            touch "$SV_DIR/.complete"
            log_event OK "results: $SV_DIR/events.tsv (per locus), $SV_DIR/sample_summary.tsv, $SV_DIR/sv_genotype.log"
        fi
        """
    )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s18c small variants between samples and references --------------------

SNV_ENGINE = Path(__file__).resolve().parent / "snv_tools.py"
SNV_CALLERS = ("C3P", "C3AF", "DP", "INT", "UNI")


def _snv_references(cfg: dict) -> List[tuple]:
    """snv_references entries: a FASTA path, or a mapping with fasta (and optional gff, name)."""
    out = []
    for e in cfg.get("snv_references") or []:
        if isinstance(e, dict):
            fa = os.path.expandvars(os.path.expanduser(str(e.get("fasta", ""))))
            gff = os.path.expandvars(os.path.expanduser(str(e.get("gff") or "")))
            nm = str(e.get("name") or "")
        else:
            fa, gff, nm = os.path.expandvars(os.path.expanduser(str(e))), "", ""
        if not nm:
            nm = Path(fa).name
            for ext in (".gz", ".fasta", ".fna", ".fa", ".fas"):
                if nm.endswith(ext):
                    nm = nm[: -len(ext)]
        out.append((nm, fa, gff))
    return out


def stage_snv(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str | None:
    """Small variants (SNPs, indels) between samples and against references (s18c).

    Written only with snv_enabled: true. Comparisons, each a mapping of one
    sample's reads to one genome: self (the sample's own consensus: residual
    errors), cross (the consensus of a near-identical sample, e.g. replicate
    or parent; groups are formed from measured identity and aligned fraction and
    the first sample of a group in samplesheet order is its reference) and
    reference (each genome in snv_references that is close enough to the
    sample). Reads: a deterministic snv_depth-fold subsample mapped with
    minimap2 -x lr:hq; caller configuration snv_caller (validated: README,
    rule R2). Every call gets strand-resolved read support, repeat and
    homopolymer flags, gene effect (the reference's GFF, the sample's annotation
    for self/cross, or the sample's annotation lifted over to an external
    reference) and whether the two assemblies also differ there.
    Outputs: snv/comparisons.tsv, snv/calls.tsv, snv/plan.tsv,
    snv/similarity.tsv, snv/<comparison>/.
    """
    if not bool(cfg.get("snv_enabled", False)):
        return None
    name = "s18c_snv"
    threads = int(cfg.get("snv_threads", 16))
    env = cfg.get("conda_env_snv", "snv")
    caller = str(cfg.get("snv_caller", SNV_DEFAULT_CALLER))
    depth = int(cfg.get("snv_depth", SNV_DEFAULT_DEPTH))
    refs = _snv_references(cfg)
    do_self = bool(cfg.get("snv_self", True))
    do_cross = bool(cfg.get("snv_cross", True))
    cross_id = float(cfg.get("snv_cross_min_identity", 0.999))
    cross_al = float(cfg.get("snv_cross_min_aligned", 0.9))
    ref_id = float(cfg.get("snv_ref_min_identity", 0.99))
    ref_al = float(cfg.get("snv_ref_min_aligned", 0.8))
    min_af = float(cfg.get("snv_clair3_min_af", 0.65))
    clair3_env = str(cfg.get("conda_env_clair3", "clair3"))
    clair3_model = os.path.expandvars(os.path.expanduser(str(cfg.get("snv_clair3_model") or "")))
    model_name = Path(str(cfg["basecall_model"])).name
    models_dir = str(cfg.get("polish_models_dir") or Path(str(cfg["basecall_model"])).parent)
    uses_clair3 = caller != "DP"
    uses_dorado = caller in ("DP", "INT", "UNI")
    problems = []
    if not SNV_ENGINE.is_file() or not POLISH_ENGINE.is_file():
        problems.append("snv_tools.py and polish_tools.py must ship next to generate_pipeline.py")
    if caller not in SNV_CALLERS:
        problems.append(f"snv_caller must be one of {', '.join(SNV_CALLERS)} (got {caller!r})")
    if not 5 <= depth <= 500:
        problems.append("snv_depth must be between 5 and 500")
    if uses_clair3 and not Path(clair3_model).is_dir():
        problems.append(f"snv_clair3_model is not a directory: {clair3_model!r}")
    missing = [fa for _, fa, _ in refs if not Path(fa).is_file()] + [g for _, _, g in refs if g and not Path(g).is_file()]
    if missing:
        problems.append("snv_references not found: " + ", ".join(missing))
    names = [n for n, _, _ in refs]
    if len(set(names)) != len(names) or set(names) & {s.sample for s in samples}:
        problems.append("snv_references names must be unique and differ from sample names")
    if not (0.9 <= cross_id <= 1 and 0.9 <= ref_id <= 1 and 0 < cross_al <= 1 and 0 < ref_al <= 1 and 0.5 <= min_af < 1):
        problems.append("snv_cross_min_identity / snv_ref_min_identity in [0.9, 1]; aligned fractions in (0, 1]; "
                        "snv_clair3_min_af in [0.5, 1)")
    if problems:
        raise SystemExit("pipeline_config.yaml (s18c_snv): " + "; ".join(problems))
    st_engine = layout.submit_dir / "_snv_tools.py"
    pt_engine = layout.submit_dir / "_polish_tools.py"
    for src, dst in ((SNV_ENGINE, st_engine), (POLISH_ENGINE, pt_engine)):
        shutil.copyfile(src, dst)
        dst.chmod(dst.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg.get("sge_h_vmem_snv", "6G"),
        threads, cfg["email"], hold_jid=hold,
    )
    body += conda_block(env) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    snv = layout.root / "snv"
    ref_lines = "".join(f"printf '%s\\t%s\\t%s\\n' \"{n}\" \"{fa}\" \"{g}\" >> \"$SNV_DIR/references.tsv\"\n"
                        for n, fa, g in refs)
    plan_flags = ("--self " if do_self else "") + ("--cross " if do_cross else "")
    body += textwrap.dedent(
        f"""
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
        set +e
        set -o pipefail
        SNV_DIR={snv}
        ST={st_engine}
        PT={pt_engine}
        DORADO={cfg['dorado_binary']}
        mkdir -p "$SNV_DIR"
        if [ -s "$SNV_DIR/comparisons.tsv" ] && [ -f "$SNV_DIR/.complete" ]; then
            log_event SKIP "results already in $SNV_DIR (delete that directory to re-run)"
            exit 0
        fi
        if ! python3 -c "import mappy, pysam, numpy" 2>/dev/null; then
            log_event WARN "the '{env}' conda env lacks mappy/pysam/numpy (needs python>=3.9, mappy, pysam, numpy, scipy, minimap2, samtools, bcftools); skipping"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            exit 0
        fi
        for tool in minimap2 samtools bcftools; do
            if ! command -v "$tool" >/dev/null 2>&1; then
                log_event WARN "$tool not found in the '{env}' env; skipping"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                exit 0
            fi
        done
        CL3=""
        if [ "{str(uses_clair3).lower()}" = true ]; then
            CL3="{clair3_env}"
            case "$CL3" in
                /*) ;;
                *) CL3=$(conda env list 2>/dev/null | awk -v n="$CL3" '$1 == n {{print $NF}}' | head -1) ;;
            esac
            if [ -z "$CL3" ] || [ ! -x "$CL3/bin/run_clair3.sh" ]; then
                log_event WARN "Clair3 (run_clair3.sh) not found in conda env '{clair3_env}' (conda_env_clair3); skipping"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                exit 0
            fi
        fi
        if [ "{str(uses_dorado).lower()}" = true ] && [ ! -x "$DORADO" ]; then
            log_event WARN "dorado not found at $DORADO (dorado_binary); skipping"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            exit 0
        fi
        : > "$SNV_DIR/samples.tsv.tmp"
        : > "$SNV_DIR/references.tsv"
        """
    )
    body += ref_lines
    for s in samples:
        cons = f"{layout.consensus}/{s.sample}.fasta"
        reads = f"{layout.fastq_filtered}/{s.sample}.fastq"
        adir = f"{layout.annotate}/{s.sample}"
        body += textwrap.dedent(
            f"""
            # ---- {s.sample} ----
            gff=""
            for g in {adir}/{s.sample}.gff3 {adir}/{s.sample}.gff; do
                if [ -s "$g" ]; then gff="$g"; break; fi
            done
            fna={adir}/{s.sample}.fna
            [ -s "$fna" ] || fna=""
            if [ ! -s {cons} ] || [ ! -s {reads} ]; then
                log_event WARN "{s.sample}: consensus or filtered reads missing; no comparisons"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "{s.sample}" {cons} {reads} "$gff" "$fna" >> "$SNV_DIR/samples.tsv.tmp"
                [ -n "$gff" ] || log_event WARN "{s.sample}: no annotation GFF; effects will read no_annotation"
            fi
            """
        )
    body += textwrap.dedent(
        f"""
        mv "$SNV_DIR/samples.tsv.tmp" "$SNV_DIR/samples.tsv"
        if [ ! -s "$SNV_DIR/samples.tsv" ]; then
            log_event WARN "no sample has a consensus and reads; nothing to compare"
            exit 0
        fi
        if ! python3 "$ST" plan --samples "$SNV_DIR/samples.tsv" --references "$SNV_DIR/references.tsv" {plan_flags}\\
                --cross-min-identity {cross_id} --cross-min-aligned {cross_al} --ref-min-identity {ref_id} \\
                --ref-min-aligned {ref_al} --out "$SNV_DIR/plan.tsv" --out-similarity "$SNV_DIR/similarity.tsv" \\
                > "$SNV_DIR/plan.log" 2>&1; then
            log_event WARN "comparison planning failed: $(tail -3 "$SNV_DIR/plan.log" | tr '\\n' ' ' | cut -c1-400)"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            exit 0
        fi
        log_event OK "planned comparisons $(tail -1 "$SNV_DIR/plan.log"); caller {caller} at {depth}x; pairwise identities in $SNV_DIR/similarity.tsv"
        while IFS= read -r cmp; do
            if [ -f "$SNV_DIR/$cmp/.done" ]; then
                log_event SKIP "$cmp (already done)"
                continue
            fi
            rm -rf "$SNV_DIR/$cmp"
            t_c=$(date +%s)
            if python3 "$ST" run --plan "$SNV_DIR/plan.tsv" --comparison "$cmp" --outdir "$SNV_DIR/$cmp" \\
                    --depth {depth} --caller {caller} --threads {threads} --af {min_af} --polish-tools "$PT" \\
                    --stats-dir "$SNV_DIR/_read_stats" --basecall-model "{model_name}" --clair3-prefix "$CL3" \\
                    --clair3-model "{clair3_model}" --dorado "$DORADO" --models-dir "{models_dir}" --device cpu \\
                    > "$SNV_DIR/$cmp.log" 2>&1; then
                log_event OK "$cmp: $(tail -1 "$SNV_DIR/$cmp.log" | cut -c1-200) in $(( $(date +%s) - t_c )) s"
            else
                log_event WARN "$cmp failed: $(tail -4 "$SNV_DIR/$cmp.log" | tr '\\n' ' ' | cut -c1-500)"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            fi
        done < <(tail -n +2 "$SNV_DIR/plan.tsv" | cut -f1)
        if python3 "$ST" collect --plan "$SNV_DIR/plan.tsv" --snv-dir "$SNV_DIR" > "$SNV_DIR/collect.log" 2>&1; then
            # awk, not read: empty TSV fields would shift bash read's variables
            awk -F'\\t' 'NR == 1 {{ for (i = 1; i <= NF; i++) h[$i] = i; next }}
                {{ m = $h["comparison"] ": " $h["n_calls"] " calls (" $h["n_snp"] " SNP, " $h["n_indel"] " indel, " $h["n_mnp_complex"] " MNP/complex; " $h["n_coding_change"] " change a protein)"
                  if ($h["kind"] != "self") m = m "; assemblies differ at " $h["assembly_diffs"] " small sites, " $h["assembly_diffs_without_call"] " of them without a call; " $h["calls_not_in_assembly_diff"] " calls where the assemblies agree"
                  if ($h["status"] != "ok") m = m " [" $h["status"] "]"
                  print m }}' "$SNV_DIR/comparisons.tsv" | while IFS= read -r line; do log_event OK "$line"; done
            [ "$STAGE_FAIL_COUNT" -eq 0 ] && touch "$SNV_DIR/.complete"
            log_event OK "results: $SNV_DIR/comparisons.tsv, $SNV_DIR/calls.tsv"
        else
            log_event WARN "collecting results failed: $(tail -3 "$SNV_DIR/collect.log" | tr '\\n' ' ' | cut -c1-300)"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
        fi
        """
    )
    # pipefail is on for the steps above; the final scan's grep finds nothing on a clean
    # run and must not turn that into the script's exit status.
    body += '\nset +o pipefail\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s18d minority variants within samples -------------------------------

MINOR_ENGINE = Path(__file__).resolve().parent / "minor_variants.py"


def stage_minor(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str | None:
    """Minority single-nucleotide variants inside each sample (s18d).

    Written only with minor_enabled: true. On each sample's s13 mapping,
    minor_variants.py counts bases per strand, learns the sample's own error
    rates per sequence context (k-mer and homopolymer classes, beta-binomial),
    tests every position on both strands (Benjamini-Hochberg over all tests),
    removes strand-biased and context-enriched candidates, masks repeats and
    contig ends, reports error-corrected fractions with 95% intervals and links
    sites carried by the same reads into groups (haplotypes). Genome versions in
    sv_variant_genomes are also genotyped at their SNP sites.
    Outputs: minor/<sample>/{sites,groups,summary,null_model}.tsv (+ known_*.tsv),
    minor/summary.tsv.
    """
    if not bool(cfg.get("minor_enabled", False)):
        return None
    name = "s18d_minor"
    env = cfg.get("conda_env_snv", "snv")
    fdr = float(cfg.get("minor_fdr", 0.01))
    min_alt = int(cfg.get("minor_min_alt_per_strand", 3))
    min_sd = int(cfg.get("minor_min_strand_depth", 20))
    min_bq = int(cfg.get("minor_min_baseq", 10))
    min_mq = int(cfg.get("minor_min_mapq", 20))
    end_mask = int(cfg.get("minor_end_mask", 500))
    max_depth = float(cfg.get("minor_max_depth", 2000))
    parallel = int(cfg.get("minor_parallel", 4))
    variants = [os.path.expandvars(os.path.expanduser(str(v))) for v in (cfg.get("sv_variant_genomes") or [])]
    problems = []
    if not MINOR_ENGINE.is_file():
        problems.append(f"{MINOR_ENGINE} is missing (it ships next to generate_pipeline.py)")
    if not 0 < fdr <= 0.1:
        problems.append("minor_fdr must be in (0, 0.1]")
    if (min_alt < 2 or min_sd < 10 or not 0 <= min_bq <= 60 or not 0 <= min_mq <= 60 or end_mask < 0
            or max_depth < 0 or parallel < 1):
        problems.append("minor_min_alt_per_strand >= 2; minor_min_strand_depth >= 10; minor_min_baseq and "
                        "minor_min_mapq in [0, 60]; minor_end_mask, minor_max_depth >= 0; minor_parallel >= 1")
    missing = [v for v in variants if not Path(v).is_file()]
    if missing:
        problems.append("sv_variant_genomes not found: " + ", ".join(missing))
    if problems:
        raise SystemExit("pipeline_config.yaml (s18d_minor): " + "; ".join(problems))
    engine = layout.submit_dir / "_minor_variants.py"
    shutil.copyfile(MINOR_ENGINE, engine)
    engine.chmod(engine.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    slots = max(1, min(parallel, len(samples)))
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg.get("sge_h_vmem_minor", "8G"),
        slots, cfg["email"], hold_jid=hold,
    )
    body += conda_block(env) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    minor = layout.root / "minor"
    ver = " ".join(f'--version "{v}"' for v in variants) if len(variants) >= 2 else ""
    opts = (f"--fdr {fdr} --min-alt-per-strand {min_alt} --min-strand-depth {min_sd} --min-baseq {min_bq} "
            f"--min-mapq {min_mq} --end-mask {end_mask} --max-depth {max_depth:g}")
    body += textwrap.dedent(
        f"""
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
        set +e
        MINOR_DIR={minor}
        ME={engine}
        mkdir -p "$MINOR_DIR"
        if ! python3 -c "import mappy, pysam, numpy, scipy" 2>/dev/null; then
            log_event WARN "the '{env}' conda env lacks mappy/pysam/numpy/scipy; skipping"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            exit 0
        fi
        """
    )
    body += textwrap.dedent(
        f"""
        # Samples run in parallel ({slots} at a time); each writes its exit code to
        # minor/<sample>.rc and is reported once all have finished.
        run_minor() {{
            local S="$1"; shift
            python3 "$ME" call --sample "$S" --outdir "$MINOR_DIR/$S" "$@" {opts} > "$MINOR_DIR/$S.log" 2>&1
            echo $? > "$MINOR_DIR/$S.rc"
        }}
        LAUNCHED=""
        """
    )
    for s in samples:
        cons = f"{layout.consensus}/{s.sample}.fasta"
        bam = f"{layout.mapping}/{s.sample}/mapped.bam"
        adir = f"{layout.annotate}/{s.sample}"
        body += textwrap.dedent(
            f"""
            # ---- {s.sample} ----
            gff=""
            for g in {adir}/{s.sample}.gff3 {adir}/{s.sample}.gff; do
                if [ -s "$g" ]; then gff="$g"; break; fi
            done
            fna={adir}/{s.sample}.fna
            [ -s "$fna" ] || fna=""
            if [ -f "$MINOR_DIR/{s.sample}/.done" ]; then
                log_event SKIP "{s.sample} (already done; delete $MINOR_DIR/{s.sample} to re-run)"
            elif [ ! -s {cons} ] || [ ! -s {bam} ] || [ ! -s {bam}.bai ]; then
                log_event WARN "{s.sample}: consensus or indexed s13 mapping (mapped.bam + .bai) missing; not assessed"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                rm -rf "$MINOR_DIR/{s.sample}" "$MINOR_DIR/{s.sample}.rc"
                while [ "$(jobs -rp | wc -l)" -ge {slots} ]; do sleep 10; done
                run_minor "{s.sample}" --bam {bam} --ref {cons} --gff "$gff" --gff-fasta "$fna" {ver} &
                LAUNCHED="$LAUNCHED {s.sample}"
                log_event OK "{s.sample}: started"
            fi
            """
        )
    body += textwrap.dedent(
        f"""
        wait
        for S in $LAUNCHED; do
            rc=$(cat "$MINOR_DIR/$S.rc" 2>/dev/null || echo missing)
            if [ "$rc" = 0 ] && [ -s "$MINOR_DIR/$S/summary.tsv" ]; then
                touch "$MINOR_DIR/$S/.done"
                log_event OK "$S: $(awk -F'\\t' 'NR == 1 {{ for (i = 1; i <= NF; i++) h[$i] = i; next }} {{ printf "%s called and %s linked minority sites in %s groups; %s bp assessed at median depth %s (read fraction used %s), median detection floor %s; %s s", $h["n_called"], $h["n_linked"], $h["n_groups"], $h["positions_assessed"], $h["median_depth"], $h["read_fraction_used"], $h["floor_median"], $h["seconds"] }}' "$MINOR_DIR/$S/summary.tsv")"
            else
                log_event WARN "$S: minority-variant engine failed (rc=$rc): $(tail -3 "$MINOR_DIR/$S.log" | tr '\\n' ' ' | cut -c1-400)"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            fi
            rm -f "$MINOR_DIR/$S.rc"
        done
        """
    )
    body += textwrap.dedent(
        f"""
        # run-level summary: one row per sample
        first=1
        : > "$MINOR_DIR/summary.tsv.tmp"
        for f in {" ".join(f'"$MINOR_DIR/{s.sample}/summary.tsv"' for s in samples)}; do
            [ -s "$f" ] || continue
            if [ $first -eq 1 ]; then cat "$f"; first=0; else tail -n +2 "$f"; fi >> "$MINOR_DIR/summary.tsv.tmp"
        done
        mv "$MINOR_DIR/summary.tsv.tmp" "$MINOR_DIR/summary.tsv"
        log_event OK "results: $MINOR_DIR/<sample>/sites.tsv, groups.tsv; run summary $MINOR_DIR/summary.tsv"
        """
    )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s18e phase variation (PhaVa) ---------------------------------------

PHAVA_ENGINE = Path(__file__).resolve().parent / "phava_tools.py"


def stage_phava(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str | None:
    """Phase variation: invertible DNA segments with PhaVa (s18e).

    Written only with phava_enabled: true. For each sample, phava_tools.py turns
    the annotation into the GFF that PhaVa reads (region lines, consensus contig
    names, transl_table, no FASTA section) and subsamples the filtered reads by
    read name to phava_max_depth; PhaVa 0.2.3 (variation_wf --reportAll) finds
    inverted repeats in the consensus with einverted and counts the reads that
    support each orientation of every candidate inverton; phava_tools.py
    genotype re-counts every candidate with an inverted read against its own two
    orientations only (PhaVa's joint mapping loses reference-orientation reads
    where candidates are close, inflating the inverted fraction), using the s13
    mapping to find the reads; phava_tools.py annotate adds the reverse fraction
    with a 95% interval, the call (at least phava_min_reverse_reads reads in the
    inverted orientation, making up at least phava_min_reverse_fraction of the
    counted reads) and the gene context.
    Outputs: phava/<sample>/{invertons,summary,genotype,phava_ratio,phava_IRs}.tsv (+ the
    GFF PhaVa used, phava_genes.gff), phava/{summary,invertons_called,invertons_all}.tsv.
    """
    if not bool(cfg.get("phava_enabled", False)):
        return None
    name = "s18e_phava"
    env = cfg.get("conda_env_phava", "phava")
    max_depth = float(cfg.get("phava_max_depth", 1000))
    min_rev = int(cfg.get("phava_min_reverse_reads", 3))
    min_frac = float(cfg.get("phava_min_reverse_fraction", 0.03))
    threads = int(cfg.get("phava_threads", 4))
    parallel = int(cfg.get("phava_parallel", 4))
    problems = []
    if not PHAVA_ENGINE.is_file():
        problems.append(f"{PHAVA_ENGINE} is missing (it ships next to generate_pipeline.py)")
    if max_depth < 0 or min_rev < 1 or not 0 <= min_frac < 1 or threads < 1 or parallel < 1:
        problems.append("phava_max_depth >= 0; phava_min_reverse_reads >= 1; phava_min_reverse_fraction "
                        "in [0, 1); phava_threads and phava_parallel >= 1")
    if problems:
        raise SystemExit("pipeline_config.yaml (s18e_phava): " + "; ".join(problems))
    engine = layout.submit_dir / "_phava_tools.py"
    shutil.copyfile(PHAVA_ENGINE, engine)
    engine.chmod(engine.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    par = max(1, min(parallel, len(samples)))
    body = sge_header(
        name, layout.stdout_dir, cfg["sge_h_rt_cpu"], cfg.get("sge_h_vmem_phava", "4G"),
        par * threads, cfg["email"], hold_jid=hold,
    )
    body += conda_block(env) + "\n\n"
    body += log_setup(name, layout) + "\n\n"
    phava = layout.root / "phava"
    opts = f"--min-reverse-reads {min_rev} --min-reverse-fraction {min_frac:g} --promoter-window 300"
    body += textwrap.dedent(
        f"""
        # numpy/OpenBLAS (imported by PhaVa) reserves buffers for every core of the node,
        # which breaks the h_vmem limit on large nodes.
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
        set +e
        PHAVA_DIR={phava}
        PT={engine}
        mkdir -p "$PHAVA_DIR"
        for t in phava einverted minimap2; do
            if ! command -v "$t" > /dev/null 2>&1; then
                log_event WARN "the '{env}' conda env lacks $t (phava 0.2.3 from bioconda brings all three); skipping"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                exit 0
            fi
        done
        PHAVA_VERSION=$(phava -h 2>&1 | grep -o 'PhaVa v[0-9.]*' | head -n 1 | cut -c8-)
        """
    )
    body += textwrap.dedent(
        f"""
        # Samples run in parallel ({par} at a time, {threads} threads each); each writes its
        # exit code to phava/<sample>.rc and is reported once all have finished. PhaVa
        # names its outputs after the read and genome files, so both are named <sample>.
        run_phava() {{
            local S="$1" CONS="$2" FASTQ="$3" GFF="$4" FNA="$5" BAM="$6"
            local D="$PHAVA_DIR/$S"
            local W="$PHAVA_DIR/$S/work"
            (
                set -e
                mkdir -p "$W"
                cp "$CONS" "$W/$S.fasta"
                GOPT=""
                if [ -n "$GFF" ]; then
                    python3 "$PT" gff --gff "$GFF" --fna "$FNA" --consensus "$W/$S.fasta" \\
                        --out "$D/phava_genes.gff" --json "$D/phava_genes.json"
                    GOPT="--genes $D/phava_genes.gff --genesFormat gff"
                fi
                python3 "$PT" subsample --fastq "$FASTQ" --genome "$W/$S.fasta" --max-depth {max_depth:g} \\
                    --seed phava1 --out "$W/$S.fastq" --json "$D/subsample.json"
                [ -e "$W/$S.fastq" ] || ln -s "$FASTQ" "$W/$S.fastq"
                phava variation_wf -i "$W/$S.fasta" -r "$W/$S.fastq" -d "$W/phava" -t {threads} $GOPT \\
                    --reportAll > "$W/phava.stdout"
                cp "$W/phava/$S.fastq_vs_$S.fasta_ratio.tsv" "$D/phava_ratio.tsv"
                cp "$W/phava/data/IRs.tsv" "$D/phava_IRs.tsv"
                cp "$W/phava/Phava.log" "$D/phava.log"
                GT=""
                if [ -n "$BAM" ]; then
                    python3 "$PT" genotype --ratio "$D/phava_ratio.tsv" --consensus "$W/$S.fasta" --bam "$BAM" \\
                        --fastq "$W/$S.fastq" --subsample-json "$D/subsample.json" --threads {threads} \\
                        --tmp "$W" --out "$D/genotype.tsv"
                    GT="--genotype $D/genotype.tsv"
                fi
                python3 "$PT" annotate --ratio "$D/phava_ratio.tsv" --consensus "$W/$S.fasta" --gff "$GFF" \\
                    --fna "$FNA" --sample "$S" {opts} --subsample-json "$D/subsample.json" $GT \\
                    --phava-version "$PHAVA_VERSION" --out "$D/invertons.tsv" --summary "$D/summary.tsv"
                rm -rf "$W"
            ) > "$D.log" 2>&1
            echo $? > "$D.rc"
        }}
        LAUNCHED=""
        """
    )
    for s in samples:
        cons = f"{layout.consensus}/{s.sample}.fasta"
        reads = f"{layout.fastq_filtered}/{s.sample}.fastq"
        adir = f"{layout.annotate}/{s.sample}"
        bam = f"{layout.mapping}/{s.sample}/mapped.bam"
        body += textwrap.dedent(
            f"""
            # ---- {s.sample} ----
            bam={bam}
            [ -s "$bam" ] && [ -s "$bam.bai" ] || bam=""
            gff=""
            for g in {adir}/{s.sample}.gff3 {adir}/{s.sample}.gff; do
                if [ -s "$g" ]; then gff="$g"; break; fi
            done
            fna={adir}/{s.sample}.fna
            [ -s "$fna" ] || fna=""
            if [ -f "$PHAVA_DIR/{s.sample}/.done" ]; then
                log_event SKIP "{s.sample} (already done; delete $PHAVA_DIR/{s.sample} to re-run)"
            elif [ ! -s {cons} ] || [ ! -s {reads} ]; then
                log_event WARN "{s.sample}: consensus or filtered reads missing; not assessed"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            else
                [ -n "$gff" ] || log_event WARN "{s.sample}: no annotation; invertons without gene context"
                [ -n "$bam" ] || log_event WARN "{s.sample}: no indexed s13 mapping; fractions are PhaVa's joint counts (not re-counted)"
                rm -rf "$PHAVA_DIR/{s.sample}" "$PHAVA_DIR/{s.sample}.rc"
                while [ "$(jobs -rp | wc -l)" -ge {par} ]; do sleep 10; done
                run_phava "{s.sample}" {cons} {reads} "$gff" "$fna" "$bam" &
                LAUNCHED="$LAUNCHED {s.sample}"
                log_event OK "{s.sample}: started"
            fi
            """
        )
    body += textwrap.dedent(
        f"""
        wait
        for S in $LAUNCHED; do
            rc=$(cat "$PHAVA_DIR/$S.rc" 2>/dev/null || echo missing)
            if [ "$rc" = 0 ] && [ -s "$PHAVA_DIR/$S/summary.tsv" ]; then
                touch "$PHAVA_DIR/$S/.done"
                log_event OK "$S: $(awk -F'\\\\t' 'NR == 1 {{ for (i = 1; i <= NF; i++) h[$i] = i; next }} {{ printf "%s candidate invertons, %s spanned by reads, %s with inverted reads, %s called (%s intergenic, %s at a gene boundary, %s intragenic) at %sx; %s re-counted alone (PhaVa joint counts: %s called)", $h["candidates"], $h["spanned"], $h["with_reverse_reads"], $h["called"], $h["called_intergenic"], $h["called_gene_boundary"], $h["called_intragenic"], $h["depth_used"], $h["regenotyped"], $h["phava_called"] }}' "$PHAVA_DIR/$S/summary.tsv")"
            else
                log_event WARN "$S: PhaVa failed (rc=$rc): $(tail -3 "$PHAVA_DIR/$S.log" | tr '\\\\n' ' ' | cut -c1-400)"
                STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
            fi
            rm -f "$PHAVA_DIR/$S.rc"
        done
        """
    )
    body += textwrap.dedent(
        f"""
        # run-level tables: summary, called invertons (linked across samples), all candidates
        if python3 "$PT" collect --dir "$PHAVA_DIR" --samples {",".join(s.sample for s in samples)} > "$PHAVA_DIR/collect.log" 2>&1; then
            log_event OK "results: $PHAVA_DIR/<sample>/invertons.tsv; run tables $PHAVA_DIR/summary.tsv, invertons_called.tsv, invertons_all.tsv"
        else
            log_event WARN "could not combine the PhaVa tables: $(tail -2 "$PHAVA_DIR/collect.log" | tr '\\\\n' ' ' | cut -c1-300)"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
        fi
        """
    )
    body += '\nlog_warn_errors_from "$STAGE_LOG"\n'
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
                    # Bakta's JSON uses lower-case type names for some features
                    # ('cds', 'sorf', 'crispr'); they never matched the keys above,
                    # so every Bakta-annotated genome reported 0 CDS.
                    t = {{"cds": "CDS", "sorf": "sORF", "crispr": "CRISPR",
                          "ncRNA-region": "ncRNA"}}.get(t, t)
                    if t in counts:
                        counts[t] += 1
                if any(counts.values()):
                    return counts
            # Prokka (or Bakta without JSON): parse the .txt summary.
            t = sample_dir / f"{{stem}}.txt"
            if t.exists():
                # Prokka keys: CDS, rRNA, tRNA, tmRNA, misc_RNA, repeat_region;
                # Bakta's .txt uses plural keys (CDSs, tRNAs, CRISPR arrays, ...)
                alias = {{"misc_RNA": "ncRNA", "repeat_region": "CRISPR",
                          "CDSs": "CDS", "tRNAs": "tRNA", "rRNAs": "rRNA",
                          "tmRNAs": "tmRNA", "ncRNAs": "ncRNA",
                          "CRISPR arrays": "CRISPR", "sORFs": "sORF",
                          "oriCs": "oriC", "oriTs": "oriT"}}
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
import argparse, csv, hashlib, os, re, sys
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


def genbank_cds_proteins(path):
    """Yield (locus_tag, gene, product, db_xrefs, translation) per CDS feature.

    Mirrors bakta.expert.protein_sequences.parse_user_protein_sequences_genbank:
    CDS features carrying /pseudo or /pseudogene are skipped; /locus_tag,
    /product and /translation are required; /db_xref values plus
    'EC:<EC_number>' become the dbxrefs.  Plain-text parser, no Biopython.
    """
    def unquote(v):
        v = v.strip()
        if len(v) >= 2 and v[0] == '"' and v[-1] == '"':
            v = v[1:-1]
        return v.replace('""', '"')

    def finish(feat):
        if feat is None or feat[0] != "CDS":
            return None
        q = feat[1]
        if "pseudo" in q or "pseudogene" in q:
            return None
        def first(k):
            return unquote(q[k][0]) if q.get(k) else ""
        locus_tag = first("locus_tag").strip()
        product = " ".join(first("product").split())
        seq = re.sub(r"[^A-Za-z*]", "", first("translation")).upper()
        if not (locus_tag and product and seq):
            return None
        gene = first("gene").strip()
        dbx = [unquote(x).strip() for x in q.get("db_xref", [])]
        dbx = [x for x in dbx if ":" in x]
        dbx += ["EC:" + unquote(e).strip() for e in q.get("EC_number", [])]
        return locus_tag, gene, product, dbx, seq

    feat, cur, in_ft = None, None, False
    with open(path, errors="ignore") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if line.startswith("FEATURES"):
                in_ft = True
                continue
            if not in_ft:
                continue
            if line[:1] not in ("", " "):
                # ORIGIN / CONTIG / BASE COUNT / '//' close the feature table
                rec = finish(feat)
                if rec:
                    yield rec
                feat, cur, in_ft = None, None, False
                continue
            if len(line) > 5 and line[:5] == "     " and line[5] != " ":
                rec = finish(feat)
                if rec:
                    yield rec
                feat, cur = (line[5:21].strip(), {}), None
                continue
            if feat is None:
                continue
            body = line[21:] if line[:21].strip() == "" else line.strip()
            if body.startswith("/"):
                name, _, val = body[1:].partition("=")
                cur = name.strip()
                feat[1].setdefault(cur, []).append(val)
            elif cur is not None:
                sep = "" if cur == "translation" else " "
                feat[1][cur][-1] += sep + body.strip()
    rec = finish(feat)
    if rec:
        yield rec


def genbank_to_bakta_fasta(p, workdir):
    """Write Bakta short-format FASTA (>locus_tag gene~~~product~~~dbxrefs)
    from the CDS translations of a GenBank file; cached per path/size/mtime."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    st = p.stat()
    key = f"{p.resolve()}:{st.st_size}:{st.st_mtime_ns}"
    out = workdir / f"{p.stem}.{hashlib.md5(key.encode()).hexdigest()[:10]}.bakta.faa"
    if out.is_file() and out.stat().st_size > 0:
        return str(out), "genbank-converted-cached"
    n = 0
    tmp = out.with_suffix(".tmp")
    with tmp.open("w") as fo:
        for locus_tag, gene, product, dbx, seq in genbank_cds_proteins(p):
            gene = gene.replace("~~~", " ")
            product = product.replace("~~~", " ")
            dbx = ",".join(x.replace(",", ";") for x in dbx)
            fo.write(f">{locus_tag} {gene}~~~{product}~~~{dbx}\n{seq}\n")
            n += 1
    if n == 0:
        tmp.unlink(missing_ok=True)
        return "", "empty"
    tmp.replace(out)
    return str(out), f"genbank-converted({n})"


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
    [Escherichia coli]') therefore has to be rewritten. GenBank input passes
    straight through to Prokka; for Bakta its CDS translations are converted
    to the short format above (genbank_to_bakta_fasta), because Bakta 1.12
    with Biopython >= 1.87 rejects GenBank given to --proteins. Prokka also
    tolerates plain FASTA.

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

    # GenBank/EMBL: Prokka reads it directly (prokka-genbank_to_fasta_db).
    # Bakta 1.12 parses --proteins as FASTA *first* and exits if that raises;
    # Biopython >= 1.87 raises on any non-FASTA text, so a GenBank file dies
    # with "User proteins file Fasta format not valid!" before Bakta's own
    # GenBank branch runs.  Convert GenBank CDS translations for Bakta.
    if p.suffix.lower() in GENBANK_SUFFIXES:
        if tool == "bakta":
            with p.open(errors="ignore") as fh:
                head = fh.readline()
            if head.startswith("LOCUS"):
                return genbank_to_bakta_fasta(p, workdir)
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
    p.add_argument("--sv-dir", type=Path, default=None,
                   help="s18b_sv output directory (events.tsv); optional")
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


    # ---- structural heterogeneity (s18b_sv) ----
    sv_html = "<p>The structural-heterogeneity stage (s18b_sv) was not run or produced no results.</p>"
    ev_path = args.sv_dir / "events.tsv" if args.sv_dir else None
    if ev_path is not None and ev_path.exists():
        try:
            ev = pd.read_csv(ev_path, sep="\t", dtype=str).fillna("")
            if ev.empty:
                sv_html = "<p>s18b_sv ran but no candidate locus could be placed in any sample.</p>"
            else:
                shown = ev[(ev["source"] == "known") | ev["status"].isin(["mixture", "consensus_discrepancy"])].copy()
                def frac_col(r):
                    if r.get("method") == "junction":
                        return f"junction: {r.get('junction_alt_fraction_ci95', '')} ({r.get('junction_call', '')})"
                    return r.get("allele_reads_fraction_ci95", "")
                def nreads(r):
                    # reads behind the call: junction observations when the junctions decided it
                    if r.get("method") == "junction" and r.get("junction_counts"):
                        try:
                            return str(sum(int(x.split("=")[1]) for x in r["junction_counts"].split("; ") if "=" in x))
                        except ValueError:
                            pass
                    return r["spanning_reads"]
                def cut(v, n=160):
                    v = str(v)
                    return v if len(v) <= n else v[:n] + "…"
                def genes_cell(v):
                    # complete list; folded when long (the TSV has the same text)
                    items = [x for x in str(v).split("; ") if x]
                    if not items:
                        return "—"
                    full = escape("; ".join(items))
                    if len(items) <= 3:
                        return full
                    return f"<details><summary>{len(items)} genes</summary>{full}</details>"
                rows_ = []
                for _, r in shown.iterrows():
                    rows_.append({
                        "sample": escape(r["sample"]), "locus": escape(r["locus"]), "kind": escape(r["kind"]),
                        "window": escape(f"{r['contig']}:{r['start']}-{r['end']} ({r['consensus_strand']})"),
                        "segment": escape(r.get("consensus_segment", "")) or "—",
                        "status": escape(r["status"]), "method": escape(r.get("method", "spanning")),
                        "reads (spanning or junction)": escape(nreads(r)),
                        "allele fractions [95% CI]": escape(cut(frac_col(r), 260)),
                        "floor": escape(r["detection_floor"]), "flags": escape(cut(r["flags"], 80)),
                        "genes in segment": genes_cell(r["genes_in_consensus_segment"]),
                        "breakpoint genes": "<br>".join(escape(x) for x in str(r.get("breakpoint_genes", "")).split(" | ") if x) or "—",
                    })
                tbl = pd.DataFrame(rows_)
                counts = ev.groupby(["sample", "status"]).size().unstack(fill_value=0).reset_index()
                sv_html = (
                    "<p class='meta'>Each candidate structural variant is re-genotyped against explicit allele "
                    "sequences: a read counts for an allele only if one alignment spans that allele's whole window "
                    "(all alleles share one window length), or, for events longer than the reads, at the breakpoint "
                    "junctions. Fractions are shares of genome copies with Wilson 95% intervals; an allele is called "
                    "present when its lower bound exceeds the sample's detection floor (highest fake-allele fraction "
                    "at size-matched null loci, at least sv_min_fraction). <em>segment</em>: the event's variable "
                    "sequence in this sample's consensus (<code>contig:p^q</code>: the junction where a sequence the "
                    "consensus lacks joins it); <em>breakpoint genes</em>: the gene at each end of the segment "
                    "([cut]: the gene extends beyond it) or at the junction, else the nearest gene on each side. "
                    f"Full per-locus results: <code>{escape(str(ev_path))}</code>.</p>"
                    "<h3>Known genome versions and reported mixtures</h3>"
                    + (tbl.to_html(index=False, classes="metrics", border=0, na_rep="—", escape=False)
                       if not tbl.empty else "<p>No known locus and no supported mixture.</p>")
                    + "<h3>All candidates by status</h3>" + table_html(counts)
                )
        except Exception as e:
            sv_html = f"<p>SV results could not be read: {escape(str(e))}</p>"

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
        "<h2>Structural heterogeneity within samples</h2>"
        f"{sv_html}"
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

REPORT_SMALLVAR_PY = r'''#!/usr/bin/env python3
"""Adds the consensus-polishing (s12b), small-variant (s18c) and minority-variant
(s18d) sections to report.html. Called by s20_report after _report.py; only the
sections of stages that are enabled are added (their directories are passed)."""
import argparse
import json
from html import escape
from pathlib import Path

import pandas as pd


def table(df, max_rows=None):
    note = ""
    if max_rows is not None and len(df) > max_rows:
        note = f"<p class='meta'>Showing the first {max_rows} of {len(df)} rows; the full table is in the file named above.</p>"
        df = df.head(max_rows)
    return df.to_html(index=False, classes="metrics", border=0, na_rep="—") + note


def read_tsv(p):
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)


def effect_class(e):
    return "intergenic" if e.startswith("intergenic") else (e or "unannotated")


def count_text(values):
    c = {}
    for v in values:
        c[v] = c.get(v, 0) + 1
    return ", ".join(f"{k} {n}" for k, n in sorted(c.items(), key=lambda x: (-x[1], x[0])))


def snv_gene_table(calls):
    """One row per comparison and gene; an intergenic call counts for its nearest gene."""
    rows = []
    for (cmp_, lt, g, pr), x in calls.groupby(["comparison", "locus_tag", "gene", "product"], sort=False):
        rows.append({"comparison": cmp_, "locus_tag": lt, "gene": g, "product": pr, "calls": len(x),
                     "effects": count_text(effect_class(e) for e in x["effect"]),
                     "protein changes": "; ".join(v for v in x["aa_change"] if v)})
    t = pd.DataFrame(rows)
    order = {c: i for i, c in enumerate(dict.fromkeys(calls["comparison"]))}
    t["_o"] = t["comparison"].map(order)
    return t.sort_values(["_o", "calls", "locus_tag"], ascending=[True, False, True], kind="stable").drop(columns="_o")


def minor_gene_table(x):
    """One row per sample and gene; an intergenic site counts for its nearest gene."""
    rows = []
    for (smp, lt, g, pr), y in x.groupby(["sample", "locus_tag", "gene", "product"], sort=False):
        rows.append({"sample": smp, "locus_tag": lt, "gene": g, "product": pr, "sites": len(y),
                     "called": int((y["status"] == "called").sum()), "linked": int((y["status"] == "linked").sum()),
                     "effects": count_text(effect_class(e) for e in y["effect"]),
                     "max fraction (%)": f"{100 * max(float(v) for v in y['fraction']):.2f}",
                     "protein changes": "; ".join(f"{a} ({100 * float(f):.1f}%)"
                                                  for a, f in zip(y["aa_change"], y["fraction"]) if a)})
    t = pd.DataFrame(rows)
    order = {s: i for i, s in enumerate(dict.fromkeys(x["sample"]))}
    t["_o"] = t["sample"].map(order)
    return t.sort_values(["_o", "sites", "locus_tag"], ascending=[True, False, True], kind="stable").drop(columns="_o")


def group_genes_from_sites(sites):
    """genes column for groups.tsv files written before it existed (same rule as the engine)."""
    per = {}
    for _, r in sites.iterrows():
        if r["group"] == "-":
            continue
        label = (f"{r['locus_tag']}{'(' + r['gene'] + ')' if r['gene'] else ''}:{r['product']}"
                 if r["locus_tag"] else "unannotated")
        cnt = per.setdefault(r["group"], {}).setdefault(label, {})
        cls = effect_class(r["effect"])
        cnt[cls] = cnt.get(cls, 0) + 1
    return {gid: " | ".join(f"{lab} [{', '.join(f'{c} x{n}' for c, n in cls.items())}]" for lab, cls in d.items())
            for gid, d in per.items()}


def polish_section(d):
    d = Path(d)
    rows, changes = [], []
    for info in sorted(d.glob("*/info.json")):
        i = json.loads(info.read_text())
        rows.append({
            "sample": i["sample"],
            "depth (target / used)": f"{i['depth_target']}x / {i['depth_realised']}x",
            "reads used": f"{int(i['reads_used']):,}",
            "changes": i["changes"],
            "by kind": ", ".join(f"{k} {v}" for k, v in sorted(i["change_kinds"].items())) or "none",
            "dorado filter": ", ".join(f"{k} {v}" for k, v in sorted(i.get("change_filters", {}).items())) or "—",
            "length draft → polished (bp)": f"{int(i['draft_length']):,} → {int(i['polished_length']):,}",
            "installed as consensus": "yes" if i.get("installed_sha256") else "no (review mode, or see pipeline.log)",
        })
        ct = info.parent / "changes.tsv"
        if ct.exists():
            c = read_tsv(ct)
            if not c.empty:
                c.insert(0, "sample", i["sample"])
                changes.append(c)
    if not rows:
        return "<p>s12b_polish produced no results (see pipeline.log).</p>"
    html = ("<p class='meta'>Each Autocycler consensus was polished with <code>dorado polish --bacteria</code> "
            "on a deterministic read subsample (validation and the choice of depth: README, rule R1). "
            "Every change is listed with the reads that support the draft (ref) and polished (alt) sequence "
            "over a window widened to whole homopolymer runs; LowQual is dorado's own low-confidence flag. "
            f"The drafts are kept in <code>{escape(str(d))}/&lt;sample&gt;/draft.fasta</code>.</p>"
            + table(pd.DataFrame(rows)))
    if changes:
        allc = pd.concat(changes, ignore_index=True)
        html += "<h3>Changes made by polishing</h3>" + table(allc, 300)
    else:
        html += "<p>Polishing changed no base in any sample.</p>"
    return html


def snv_section(d):
    d = Path(d)
    cp = d / "comparisons.tsv"
    if not cp.exists():
        return "<p>s18c_snv produced no results (see pipeline.log).</p>"
    comps = read_tsv(cp)
    show = comps[["comparison", "kind", "sample", "reference", "depth", "n_calls", "n_snp", "n_indel",
                  "n_mnp_complex", "n_coding_change", "assembly_diffs", "assembly_diffs_without_call",
                  "calls_not_in_assembly_diff", "identity", "status"]]
    html = ("<p class='meta'>Reads of each sample mapped to its own consensus (self: residual consensus errors), "
            "to the consensus of a near-identical sample (cross) and to supplied reference genomes "
            "(reference). Caller configuration and depth were chosen by validation (README, rule R2). "
            "<em>assembly_diffs</em> counts small differences between the two assemblies themselves; "
            "a difference without a read-based call, or a call where the assemblies agree, marks a site "
            f"to inspect. Files: <code>{escape(str(cp))}</code>, <code>{escape(str(d / 'calls.tsv'))}</code>.</p>"
            + table(show))
    calls = read_tsv(d / "calls.tsv") if (d / "calls.tsv").exists() else pd.DataFrame()
    if calls.empty:
        html += "<p>No call in any comparison.</p>"
        return html
    calls = calls.assign(**{"ALT reads (fwd/rev)": calls["alt_fwd"] + "/" + calls["alt_rev"],
                            "REF reads (fwd/rev)": calls["ref_fwd"] + "/" + calls["ref_rev"]})
    cols = ["comparison", "contig", "pos", "ref", "alt", "vtype", "callers", "ALT reads (fwd/rev)",
            "REF reads (fwd/rev)", "in_repeat", "hp_len", "locus_tag", "gene", "product", "effect", "aa_change",
            "assembly_diff"]
    kind = dict(zip(comps["comparison"], comps["kind"]))
    calls["_k"] = calls["comparison"].map(lambda c: {"cross": 0, "self": 1, "reference": 2}.get(kind.get(c), 3))
    calls = calls.sort_values(["_k", "comparison"], kind="stable")
    html += ("<h3>Genes with calls</h3><p class='meta'>One row per comparison and gene (an intergenic call counts "
             "for the nearest gene, named in its effect).</p>" + table(snv_gene_table(calls), 5000))
    html += "<h3>Calls (cross-sample and self comparisons first)</h3>" + table(calls[cols], 5000)
    return html


def minor_section(d):
    d = Path(d)
    sp = d / "summary.tsv"
    if not sp.exists():
        return "<p>s18d_minor produced no results (see pipeline.log).</p>"
    summ = read_tsv(sp)
    keep = [c for c in ["sample", "positions_assessed", "positions_masked", "median_depth", "floor_p10",
                        "floor_median", "floor_p90", "n_called", "n_linked", "n_majority", "n_groups",
                        "read_fraction_used"] if c in summ.columns]
    html = ("<p class='meta'>Minority single-nucleotide variants inside each sample from the s13 mapping: "
            "per-strand counts tested against the sample's own error model per sequence context, "
            "Benjamini-Hochberg FDR across all positions, strand-bias and context filters, repeats and contig "
            "ends masked. Fractions are error-corrected with 95% intervals; the detection floor is the "
            "smallest fraction that would be called at a typical position. Sites carried by the same reads "
            "form groups (status <em>linked</em>: below the call threshold alone, but on the reads of a group). "
            "Validation: README, rule R3.</p>" + table(summ[keep]))
    sites, groups = [], []
    for s in summ["sample"]:
        f = d / s / "sites.tsv"
        if f.exists():
            x = read_tsv(f)
            if not x.empty:
                x.insert(0, "sample", s)
                sites.append(x)
        g = d / s / "groups.tsv"
        if g.exists():
            y = read_tsv(g)
            if not y.empty:
                if "genes" not in y.columns:
                    y["genes"] = y["group"].map(group_genes_from_sites(read_tsv(f)) if f.exists() else {}).fillna("")
                y.insert(0, "sample", s)
                groups.append(y)
    if sites:
        x = pd.concat(sites, ignore_index=True)
        x = x.assign(**{"fraction [95% CI]": x["fraction"] + " [" + x["ci_low"] + "–" + x["ci_high"] + "]",
                        "ALT reads (fwd/rev)": x["alt_fwd"] + "/" + x["alt_rev"],
                        "depth (fwd/rev)": x["depth_fwd"] + "/" + x["depth_rev"]})
        cols = ["sample", "contig", "pos", "ref", "alt", "status", "fraction [95% CI]", "ALT reads (fwd/rev)",
                "depth (fwd/rev)", "q", "context", "hp_len", "group", "locus_tag", "gene", "product", "effect",
                "aa_change"]
        html += ("<h3>Genes with minority variants</h3><p class='meta'>One row per sample and gene (an "
                 "intergenic site counts for the nearest gene, named in its effect).</p>" + table(minor_gene_table(x), 5000))
        html += "<h3>Minority sites</h3>" + table(x[[c for c in cols if c in x.columns]], 5000)
    else:
        html += "<p>No minority site in any sample.</p>"
    if groups:
        gt = pd.concat(groups, ignore_index=True)
        gcols = [c for c in gt.columns if c not in ("genes", "positions")] + [c for c in ("genes", "positions") if c in gt.columns]
        html += "<h3>Linked groups</h3>" + table(gt[gcols], 1000)
    return html


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--polish-dir", default="")
    ap.add_argument("--snv-dir", default="")
    ap.add_argument("--minor-dir", default="")
    a = ap.parse_args()
    html = a.report.read_text()
    parts = []
    for title, d, fn in (("Consensus polishing", a.polish_dir, polish_section),
                         ("Small variants between samples and references", a.snv_dir, snv_section),
                         ("Minority variants within samples", a.minor_dir, minor_section)):
        if not d:
            continue
        try:
            sec = fn(d)
        except Exception as e:  # a broken table must not break the whole report
            sec = f"<p>Results in {escape(str(d))} could not be read: {escape(str(e))}</p>"
        parts.append(f"<h2>{escape(title)}</h2>{sec}")
    block = "".join(parts)
    for marker in ("<h2>Plots</h2>", "</div></body></html>"):
        if marker in html:
            html = html.replace(marker, block + marker, 1)
            break
    else:
        html += block
    a.report.write_text(html)
    print(f"added {len(parts)} section(s) to {a.report}")


if __name__ == "__main__":
    main()
'''


REPORT_PHAVA_PY = r'''#!/usr/bin/env python3
"""Adds the phase-variation (s18e_phava, PhaVa) section to report.html. Called by
s20_report after _report.py (and _report_smallvar.py); written only when s18e_phava
is enabled. With --minor-dir, the s18d minority-variant groups that overlap a
candidate inverton are listed (an inverted segment shows up there as a cluster of
linked mismatches)."""
import argparse
from html import escape
from pathlib import Path

import pandas as pd


def table(df, max_rows=None):
    note = ""
    if max_rows is not None and len(df) > max_rows:
        note = f"<p class='meta'>Showing the first {max_rows} of {len(df)} rows; the full table is in the file named above.</p>"
        df = df.head(max_rows)
    return df.to_html(index=False, classes="metrics", border=0, na_rep="—") + note


def read_tsv(p):
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)


def pct_ci(r):
    if r["reverse_fraction"] in ("", "NA"):
        return "NA"
    return (f"{100 * float(r['reverse_fraction']):.1f} [{100 * float(r['ci_low']):.1f}–"
            f"{100 * float(r['ci_high']):.1f}]")


def span(r):
    return f"{r['contig']}:{int(r['start']):,}–{int(r['end']):,}"


def phava_joint(r):
    f, v = int(r["phava_forward_reads"]), int(r["phava_reverse_reads"])
    return f"{v}/{f + v}" + (f" ({100 * v / (f + v):.1f}%)" if f + v else "")


def flank(r, side):
    g = r[f"{side}_gene"]
    if not g or g == "NA":
        return g
    return f"{g} ({r[side + '_strand']}, {r[side + '_distance']} bp)"


def gene_rows(called):
    """One row per called inverton and gene it may act on: genes overlapping the
    invertible segment, promoter targets, and for other intergenic invertons the
    nearest gene on each side."""
    rows = []
    for _, r in called.iterrows():
        rel = []
        if r["class"] != "NA":
            for g in r["genes_in_segment"].split("; "):
                if g:
                    rel.append((g, "in segment"))
            targets = [t for t in r["promoter_targets"].split("; ") if t]
            for t in targets:
                label, _, rest = t.rpartition(" (")
                rel.append((label, "promoter target (" + rest))
            if r["class"] == "intergenic" and not targets:
                for side in ("left", "right"):
                    if r[f"{side}_gene"]:
                        rel.append((r[f"{side}_gene"], f"{side} flank ({r[side + '_strand']}, "
                                    f"{r[side + '_distance']} bp)"))
        for g, how in rel:
            rows.append({"sample": r["sample"], "gene": g, "relation": how, "inverton": span(r),
                         "class": r["class"], "inverted % [95% CI]": pct_ci(r),
                         "in other samples": r.get("in_other_samples", "")})
    return pd.DataFrame(rows, columns=["sample", "gene", "relation", "inverton", "class", "inverted % [95% CI]",
                                       "in other samples"])


def minority_overlaps(allc, minor_dir, samples):
    """Every s18d group with the candidate inverton (same sample) that holds most of its sites
    (ties: called first, then the higher inverted fraction, then the leftmost); reads carrying
    an inverted segment differ from the consensus at many linked positions inside it."""
    rank = {"called": 0, "below_filter": 1, "no_reverse_reads": 2, "not_spanned": 3}
    rows = []
    for s in samples:
        g = Path(minor_dir) / s / "groups.tsv"
        if not g.exists():
            continue
        cand = allc[allc["sample"] == s]
        for _, x in read_tsv(g).iterrows():
            pos = [int(p) for p in x["positions"].split(",") if p]
            gs, ge = int(x["start"]), int(x["end"])
            best = None
            for _, r in cand[(cand["contig"] == x["contig"]) & (cand["start"].astype(int) <= ge)
                             & (cand["end"].astype(int) >= gs)].iterrows():
                k = sum(int(r["start"]) <= p <= int(r["end"]) for p in pos)
                fr = float(r["reverse_fraction"]) if r["reverse_fraction"] != "NA" else -1.0
                key = (-k, rank[r["status"]], -fr, int(r["start"]))
                if best is None or key < best[0]:
                    best = (key, k, r)
            row = {"sample": s, "minority group": x["group"], "group span": f"{x['contig']}:{gs:,}–{ge:,}",
                   "sites": x["n_sites"], "group fraction %": f"{100 * float(x['fraction']):.1f}"}
            if best is None:
                row.update({"sites inside inverton": f"0/{len(pos)}", "inverton": "", "status": "",
                            "inverted % [95% CI]": "", "PhaVa joint inverted/counted": ""})
            else:
                r = best[2]
                row.update({"sites inside inverton": f"{best[1]}/{len(pos)}", "inverton": span(r), "status": r["status"],
                            "inverted % [95% CI]": pct_ci(r), "PhaVa joint inverted/counted": phava_joint(r)})
            row["genes"] = x["genes"] if "genes" in x else ""
            rows.append(row)
    return pd.DataFrame(rows)


def phava_section(d, minor_dir=""):
    d = Path(d)
    sp = d / "summary.tsv"
    if not sp.exists():
        return "<p>s18e_phava produced no results (see pipeline.log).</p>"
    summ = read_tsv(sp)
    if summ.empty:
        return "<p>s18e_phava produced no results (see pipeline.log).</p>"
    s0 = summ.iloc[0]
    html = ("<p class='meta'>Invertible DNA segments (invertons) between inverted repeats, found in each consensus by "
            f"PhaVa {escape(s0['phava_version'])} (einverted; segments of at least 30 bp) and quantified from reads "
            "that span the whole inverton in one orientation (MAPQ above 2, at most 15% mismatches). PhaVa maps the "
            "reads to all candidates at once and then loses reference-orientation reads where candidates are close, "
            "which inflates the inverted fraction; every candidate with an inverted read is therefore re-counted "
            "against its own two orientations only (same reads and rules; PhaVa's joint counts are shown alongside). "
            "Reads were subsampled by read name where deeper than the cap (column depth used). An inverton is <em>called</em> "
            f"when at least {escape(s0['min_reverse_reads'])} reads carry the inverted orientation and they are at "
            f"least {100 * float(s0['min_reverse_fraction']):g}% of the counted reads; fractions with Wilson 95% "
            "intervals. Class is the segment against the CDS (intragenic, gene boundary, intergenic); promoter "
            "targets are flanking genes whose 5′ end faces an intergenic inverton within 300 bp. Full tables: "
            f"{escape(str(d))}/&lt;sample&gt;/invertons.tsv (all candidates), invertons_called.tsv. "
            "Validation: README, rule R4.</p>")
    st = summ.rename(columns={"candidates": "candidate invertons", "spanned": "spanned by reads",
                              "with_reverse_reads": "with inverted reads", "regenotyped": "re-counted",
                              "phava_called": "PhaVa joint counts would call", "called_intergenic": "intergenic",
                              "called_gene_boundary": "gene boundary", "called_intragenic": "intragenic",
                              "called_with_promoter_target": "with promoter target",
                              "median_counted_reads": "median counted reads", "depth_used": "depth used (x)"})
    html += table(st[["sample", "candidate invertons", "spanned by reads", "with inverted reads", "re-counted",
                      "PhaVa joint counts would call", "called",
                      "intergenic", "gene boundary", "intragenic", "with promoter target", "median counted reads",
                      "depth used (x)"]])
    cp = d / "invertons_called.tsv"
    called = read_tsv(cp) if cp.exists() else pd.DataFrame()
    if called.empty:
        html += "<p>No inverton is called in any sample.</p>"
    else:
        ct = pd.DataFrame({"sample": called["sample"], "inverton": called.apply(span, axis=1),
                           "segment (bp)": called["segment_bp"], "repeat (bp)": called["ir_bp"],
                           "inverted/counted reads": called["reverse_reads"] + "/" + called["counted_reads"],
                           "inverted % [95% CI]": called.apply(pct_ci, axis=1),
                           "PhaVa joint inverted/counted": called.apply(phava_joint, axis=1), "class": called["class"],
                           "genes in segment": called["genes_in_segment"], "breakpoint genes": called["breakpoint_genes"],
                           "promoter targets": called["promoter_targets"],
                           "left gene": called.apply(lambda r: flank(r, "left"), axis=1),
                           "right gene": called.apply(lambda r: flank(r, "right"), axis=1),
                           "in other samples": called["in_other_samples"]})
        html += ("<h3>Genes at called invertons</h3><p class='meta'>One row per called inverton and gene: genes "
                 "overlapping the invertible segment, promoter targets, and for other intergenic invertons the nearest "
                 "gene on each side.</p>" + table(gene_rows(called), 5000))
        html += "<h3>Called invertons</h3>" + table(ct, 5000)
    ap = d / "invertons_all.tsv"
    allc = read_tsv(ap) if ap.exists() else pd.DataFrame()
    below = allc[allc["status"] == "below_filter"] if not allc.empty else allc
    if not below.empty:
        bt = pd.DataFrame({"sample": below["sample"], "inverton": below.apply(span, axis=1),
                           "segment (bp)": below["segment_bp"],
                           "inverted/counted reads": below["reverse_reads"] + "/" + below["counted_reads"],
                           "inverted % [95% CI]": below.apply(pct_ci, axis=1),
                           "PhaVa joint inverted/counted": below.apply(phava_joint, axis=1), "class": below["class"],
                           "genes in segment": below["genes_in_segment"],
                           "left gene": below.apply(lambda r: flank(r, "left"), axis=1),
                           "right gene": below.apply(lambda r: flank(r, "right"), axis=1)})
        html += ("<h3>Inverted reads below the call threshold</h3><p class='meta'>Candidates with inverted reads "
                 "that do not reach the call thresholds (fewer reads or a lower fraction).</p>" + table(bt, 1000))
    if minor_dir and not allc.empty:
        ov = minority_overlaps(allc, minor_dir, list(summ["sample"]))
        if not ov.empty:
            k = ov["sites inside inverton"].map(lambda v: [int(t) for t in v.split("/")])
            half = k.map(lambda v: v[1] > 0 and 2 * v[0] >= v[1])
            n_half = int(half.sum())
            n_called = int((half & (ov["status"] == "called")).sum())
            html += ("<h3>Minority-variant groups and candidate invertons</h3><p class='meta'>Every s18d minority group "
                     "with the candidate inverton holding most of its sites. "
                     f"{n_half} of {len(ov)} groups have at least half of their sites inside a candidate inverton "
                     f"({n_called} inside a called one). Reads carrying an inverted segment differ from the consensus "
                     "at many linked positions inside it, so such a group may be an inversion rather than point "
                     "mutations. A group inside a candidate without inverted reads is not explained by that "
                     "candidate.</p>" + table(ov, 1000))
    return html


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--phava-dir", required=True)
    ap.add_argument("--minor-dir", default="")
    a = ap.parse_args()
    html = a.report.read_text()
    try:
        sec = phava_section(a.phava_dir, a.minor_dir)
    except Exception as e:  # a broken table must not break the whole report
        sec = f"<p>Results in {escape(str(a.phava_dir))} could not be read: {escape(str(e))}</p>"
    block = f"<h2>Phase variation</h2>{sec}"
    for marker in ("<h2>Plots</h2>", "</div></body></html>"):
        if marker in html:
            html = html.replace(marker, block + marker, 1)
            break
    else:
        html += block
    a.report.write_text(html)
    print(f"added the phase variation section to {a.report}")


if __name__ == "__main__":
    main()
'''


def stage_report(cfg: dict, layout: Layout, hold: str, smallvar: dict | None = None,
                 phava: dict | None = None) -> str:
    """smallvar: {"polish": dir, "snv": dir, "minor": dir} for the enabled s12b/s18c/s18d
    stages; their sections are added to report.html by _report_smallvar.py.
    phava: {"phava": dir, "minor": dir or ""} when s18e is enabled; its section is
    added by _report_phava.py."""
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
            --plots-dir {layout.reports}/plots \\
            --sv-dir {layout.sv}
        rc=$?
        set -e
        if [ $rc -ne 0 ]; then
            log_event WARN "report generator rc=$rc"
            STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
        else
            log_event OK "report written to {layout.reports}/report.html"
        fi
        {{SMALLVAR_BLOCK}}
        log_warn_errors_from "$STAGE_LOG"

        log_event SUMMARY "pipeline finished — log at $PIPELINE_LOG, summary at {layout.root}/pipeline_summary.txt, report at {layout.reports}/report.html"
        # Stage summary from pipeline_status.sh (latest attempt of each stage),
        # rewritten on every report run. It is not appended to pipeline.log: copying
        # earlier START/DONE/FAIL events back into the log re-reported them (old
        # FAILs included) on every re-run.
        PIPELINE_STATUS_SKIP="$JOB_NAME" bash "{layout.root}/pipeline_status.sh" > "{layout.root}/pipeline_summary.txt" 2>&1 || true
        cat "{layout.root}/pipeline_summary.txt"
        """
    ).strip() + "\n"
    block = ""
    if smallvar:
        sv_py = layout.submit_dir / "_report_smallvar.py"
        sv_py.write_text(REPORT_SMALLVAR_PY)
        sv_py.chmod(sv_py.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        args = " ".join(f"--{k}-dir {v}" for k, v in smallvar.items())
        block = textwrap.dedent(
            f"""
            if [ $rc -eq 0 ]; then
                if python3 {sv_py} --report {layout.reports}/report.html {args}; then
                    log_event OK "added {', '.join(smallvar)} sections to the report"
                else
                    log_event WARN "could not add the {', '.join(smallvar)} sections to the report"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                fi
            fi
            """
        ).strip()
        block = "\n" + block + "\n"
    if phava:
        ph_py = layout.submit_dir / "_report_phava.py"
        ph_py.write_text(REPORT_PHAVA_PY)
        ph_py.chmod(ph_py.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        args = f"--phava-dir {phava['phava']}" + (f" --minor-dir {phava['minor']}" if phava.get("minor") else "")
        pblock = textwrap.dedent(
            f"""
            if [ $rc -eq 0 ]; then
                if python3 {ph_py} --report {layout.reports}/report.html {args}; then
                    log_event OK "added the phase variation section to the report"
                else
                    log_event WARN "could not add the phase variation section to the report"
                    STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                fi
            fi
            """
        ).strip()
        block = (block or "\n") + pblock + "\n"
    # Without the new stages the script is byte-identical to earlier versions.
    body = body.replace("{SMALLVAR_BLOCK}", block)
    write_script(layout.submit_dir / f"{name}.sh", body)
    return name


# ---- s21 MIJAMP preprocess ----------------------------------------------

def stage_mijamp(cfg: dict, layout: Layout, samples: List[Sample], hold: str) -> str:
    """Run MIJAMP's preprocess script per sample using the methylation-tagged
    BAM from s02b_demux_bam and the corresponding Bakta-annotated .fna from
    s16_annotate. MIJAMP is expected to be git-cloned at cfg['mijamp_dir'].

    Each BAM is first downsampled to cfg['mijamp_target_depth'] (default
    100x; 0 = all reads), and cfg['mijamp_filter_threshold'] (default 0.9)
    is passed to modkit. MIJAMP's motif step keeps only sites that are
    >= 90% modified with Ndiff == 0 (no reads showing another base) and
    needs more than 100 of them. Those filters were set for ~100x data: at
    several hundred x nearly every site fails on Ndiff, and with MIJAMP's
    default threshold of 0.6 many fully methylated 6mA sites stay below 90%.
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
    target_depth = int(float(cfg.get("mijamp_target_depth", 100) or 0))
    filter_thr = float(cfg.get("mijamp_filter_threshold", 0.9))
    if target_depth < 0 or not 0 < filter_thr <= 1:
        raise SystemExit(
            "pipeline_config.yaml: mijamp_target_depth must be >= 0 and "
            f"mijamp_filter_threshold in (0, 1] (got {target_depth}, {filter_thr})")
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
                    # MIJAMP motif's site filters are absolute read counts tuned on ~100x
                    # data (a site needs Ndiff == 0: no reads showing another base). At
                    # several hundred x and above nearly every site has a few error reads,
                    # so motif finds no motifs. Downsample to mijamp_target_depth first.
                    mij_bam="$bam_file"
                    sub_bam="{layout.mijamp}/{s.sample}.downsampled.tmp.bam"
                    rm -f "$sub_bam"
                    if [ {target_depth} -gt 0 ]; then
                        glen=$(awk '!/^>/ {{ n += length($0) }} END {{ print n + 0 }}' "$genome")
                        bases=$(samtools view -@ {threads} "$bam_file" | awk -F'\\t' '{{ n += length($10) }} END {{ print n + 0 }}')
                        depth=$(awk -v b="$bases" -v g="$glen" 'BEGIN {{ printf "%.0f", (g > 0 ? b / g : 0) }}')
                        sarg=$(awk -v b="$bases" -v g="$glen" -v t={target_depth} 'BEGIN {{ if (b <= 0 || g <= 0 || t * g >= b) exit; p = int(1000000 * t * g / b + 0.5); if (p < 1) p = 1; printf "42.%06d", p }}')
                        if [ "$bases" -le 0 ]; then
                            log_event WARN "{s.sample}: could not count bases in $bam_file — not downsampling"
                        elif [ -n "$sarg" ]; then
                            log_event OK "{s.sample}: read depth ${{depth}}x, downsampling to ~{target_depth}x for MIJAMP (samtools view -s $sarg)"
                            if samtools view -@ {threads} -b -s "$sarg" -o "$sub_bam" "$bam_file"; then
                                mij_bam="$sub_bam"
                            else
                                log_event WARN "{s.sample}: downsampling failed — using the full-depth BAM"
                                rm -f "$sub_bam"
                            fi
                        else
                            log_event OK "{s.sample}: read depth ${{depth}}x <= {target_depth}x, not downsampling"
                        fi
                    fi
                    log_event OK "{s.sample}: MIJAMP preprocess (bam=$mij_bam, genome=$genome [$gsrc], modkit filter threshold {filter_thr})"
                    mkdir -p {out_dir}
                    set +e
                    "$MIJAMP_DIR/scripts/preprocess" \\
                        -b "$mij_bam" \\
                        -g "$genome" \\
                        -t {threads} \\
                        -y {filter_thr} \\
                        -o {out_dir}
                    rc=$?
                    set -e
                    rm -f "$sub_bam"
                    if [ $rc -ne 0 ]; then
                        log_event WARN "{s.sample}: MIJAMP rc=$rc"
                        STAGE_FAIL_COUNT=$((STAGE_FAIL_COUNT+1))
                    else
                        log_event OK "{s.sample}: MIJAMP preprocessing done (motif discovery is not run by this pipeline)"
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
    in-flight, stalled, or done. Considers only the latest attempt of each
    stage (its events from its last START onwards) so a re-run does not
    re-report events from earlier attempts, and reports orphans (started
    but never finished — could mean still running OR killed by SGE;
    SIGKILL bypasses the EXIT trap). Skips the "pipeline summary" blocks
    older s20_report versions appended to pipeline.log. s20_report writes
    its output to pipeline_summary.txt at the end of each run."""
    body = textwrap.dedent(
        f"""
        #!/bin/bash
        # pipeline_status.sh — quick check on what each stage is doing.
        #
        # Usage:  bash {layout.root}/pipeline_status.sh
        #
        # Reads pipeline.log and reports, for the latest attempt of each stage
        # (its events from its last [START] onwards):
        #   - stages with [START] but no [DONE]     (running, failed or killed)
        #   - stages whose [DONE] reports per-sample failures
        #   - stage-level [FAIL] events
        #   - per-level event counts
        # Earlier attempts of re-run stages stay in pipeline.log as history but are
        # not reported again.
        LOG="{layout.root}/pipeline.log"
        if [ ! -f "$LOG" ]; then
            echo "No pipeline.log at $LOG"
            exit 1
        fi

        # Older s20_report versions appended "pipeline summary" blocks (copies of
        # earlier START/DONE/FAIL events) to pipeline.log; both readers skip them.
        latest_events() {{
            awk '
                FNR == 1 {{ skip = 0 }}
                /^=+ pipeline summary =+$/ {{ skip = 1; next }}
                skip {{ if ($0 ~ /^=+$/) skip = 0; next }}
                NR == FNR {{ if ($2 == "[START]") last[$3] = FNR; next }}
                !($3 in last) || FNR >= last[$3]
            ' "$LOG" "$LOG"
        }}
        earlier_attempts() {{
            awk '
                /^=+ pipeline summary =+$/ {{ skip = 1; next }}
                skip {{ if ($0 ~ /^=+$/) skip = 0; next }}
                $2 == "[START]" && $3 != "submit_all" {{ n++; seen[$3] = 1 }}
                END {{ k = 0; for (s in seen) k++; print n - k }}
            ' "$LOG"
        }}
        EVENTS=$(latest_events)

        echo "=== event counts (latest attempt of each stage) ==="
        printf '%s\\n' "$EVENTS" | grep -aoE '\\[(START|DONE|FAIL|OK|SKIP|WARN|INFO|SCAN|SUMMARY)\\]' \\
            | sort | uniq -c | sort -rn

        echo ""
        echo "=== stages started but not yet DONE ==="
        # submit_all logs only [START]; PIPELINE_STATUS_SKIP names further stages to
        # leave out (s20_report passes itself when it writes pipeline_summary.txt).
        pending=$(printf '%s\\n' "$EVENTS" | awk -v skip=" submit_all ${{PIPELINE_STATUS_SKIP:-}} " '
            $2 == "[START]" {{ p[$3] = 1 }}
            $2 == "[DONE]"  {{ delete p[$3] }}
            END {{ for (s in p) if (index(skip, " " s " ") == 0) print "  " s }}
        ' | sort)
        echo "${{pending:-  (none)}}"

        echo ""
        echo "=== stages with per-sample failures ==="
        printf '%s\\n' "$EVENTS" | grep -aE '\\[DONE\\].*per-sample failures' || echo "  (none)"

        echo ""
        echo "=== stage-level FAILs ==="
        printf '%s\\n' "$EVENTS" | grep -aE '\\[FAIL\\]' || echo "  (none)"

        echo ""
        echo "($(earlier_attempts) earlier stage attempts not shown; full history in $LOG)"
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
        layout.mijamp, layout.sv,
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

    split_jobs = basecall_split_jobs(cfg)
    if split_jobs > 1:
        # Split basecalling: one GPU array task per chunk of pod5 files, then a
        # CPU job checks the chunk BAMs and concatenates them into basecalls.bam.
        chunks, n_kept = plan_basecall_chunks(
            layout.pod5_dir, split_jobs, layout.submit_dir / "s01_basecall_chunks.tsv")
        n01 = stage_basecall_split(cfg, layout, chunks)
        n01b = stage_merge_basecalls(cfg, layout, chunks, hold=n01)
        report_basecall_plan(cfg, chunks, split_jobs, n_kept)
    else:
        n01 = stage_basecall(cfg, layout)
        n01b = None
    n_bam = n01b or n01
    n02 = stage_demux(cfg, layout, hold=n_bam)
    # Parallel branch: dorado demux to BAM (preserves MM/ML methylation tags).
    # Runs in parallel with s02_demux off the same basecalls.bam, then
    # feeds s21_mijamp at the end.
    n02b = stage_demux_bam(cfg, layout, hold=n_bam)
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
    # Optional polishing of each consensus before anything maps to or annotates it.
    n12b = stage_polish(cfg, layout, samples, hold=n12)
    # Mapping → meta-assembly → binning happen BEFORE Bakta now, so Bakta
    # can annotate both the autocycler consensuses AND the metabat2 bins in
    # a single stage.
    n13 = stage_map_reads(cfg, layout, samples, hold=n12b or n12)
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
    # Structural heterogeneity needs the s13 read mappings and the annotation
    # (for genes in variable segments); s19/s20 wait for it so the report has it.
    n18b = stage_sv(cfg, layout, samples, hold=f"{n13},{n18}")
    # Optional small-variant stages: between samples / references, and minority
    # variants within samples. Both need the s13 mapping and the annotation.
    n18c = stage_snv(cfg, layout, samples, hold=f"{n13},{n18}")
    n18d = stage_minor(cfg, layout, samples, hold=f"{n13},{n18}")
    # Phase variation (PhaVa) maps the filtered reads itself; it needs the
    # consensus, the annotation (gene context) and the s13 mapping (re-count).
    n18e = stage_phava(cfg, layout, samples, hold=f"{n13},{n18}")
    # Aggregate reads annotation features + CheckM2 quality + GTDB-Tk taxonomy.
    n19 = stage_aggregate(cfg, layout, samples, hold=",".join([n17, n18, n18b] + [x for x in (n18c, n18d, n18e) if x]))
    smallvar = {k: str(layout.root / d) for k, d, n in (("polish", "polish", n12b), ("snv", "snv", n18c),
                                                       ("minor", "minor", n18d)) if n}
    phava = {"phava": str(layout.root / "phava"), "minor": str(layout.root / "minor") if n18d else ""} if n18e else None
    n20 = stage_report(cfg, layout, hold=n19, smallvar=smallvar or None, phava=phava)
    # MIJAMP needs both the demux BAMs (s02b) and the Bakta .fna (s16, which
    # is upstream of s20). Hold on both so SGE waits for the slower of the two.
    n21 = stage_mijamp(cfg, layout, samples, hold=f"{n20},{n02b}")

    ordered = [
        [n01],
    ] + ([[n01b]] if n01b else []) + [
        [n02], [n02b], [n03], [n04], [n05],
        n06_list,
        [n07], [n08], [n09], [n10], [n11], [n12],
    ] + ([[n12b]] if n12b else []) + [
        [n13], [n14], [n15], [n16], [n17], [n18], [n18b],
    ] + [[n] for n in (n18c, n18d, n18e) if n] + [
        [n19], [n20], [n21],
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
