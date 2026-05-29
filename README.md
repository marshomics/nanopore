# Nanopore assembly + annotation pipeline

End-to-end SGE pipeline that takes raw POD5 files and a sample sheet, runs
basecalling → demux → QC → Autocycler 8-assembler consensus → Bakta →
read-mapping → unmapped-read meta-assembly (Flye --meta) → MetaBAT2
binning → CheckM2 → GTDB-Tk, and produces an HTML report with
publication-ready plots and per-sample summary tables.

## Files

- `generate_pipeline.py` — emits all SGE submit scripts for a given run
- `pipeline_config.yaml` — paths, models, conda envs, SGE resources, filter thresholds
- `samplesheet_template.tsv` — fill in one row per barcode

## One-shot install (all environments)

Every stage runs inside a conda env named in `pipeline_config.yaml`. The
block below creates them all with mamba. You can run it on a fresh
cluster account in one go; each `mamba create` is independent and idempotent
in the sense that re-running with `--yes` upgrades in place. Adjust channel
priority or pin versions to match what your cluster admins want.

```bash
# Pick a location for non-conda software (dorado binary, MIJAMP clone, dbs)
SOFTWARE_DIR=/ebio/abt3_scratch/$USER/software
DB_DIR=/ebio/abt3_scratch/$USER/dbs
mkdir -p "$SOFTWARE_DIR" "$DB_DIR"

# --- Channels ---
# (use --override-channels in calls below to keep the env clean)
CHANNELS="-c bioconda -c conda-forge"

# === 1. nanopore — dorado (external), chopper, flye ===
# Dorado is distributed as a Linux binary by ONT, not on conda.
# Pick the version that matches the basecall model you're using.
DORADO_VER=1.3.1
cd "$SOFTWARE_DIR"
wget "https://cdn.oxfordnanoportal.com/software/analysis/dorado-${DORADO_VER}-linux-x64.tar.gz"
tar -xzf "dorado-${DORADO_VER}-linux-x64.tar.gz"
rm "dorado-${DORADO_VER}-linux-x64.tar.gz"
# dorado now lives at $SOFTWARE_DIR/dorado-${DORADO_VER}-linux-x64/bin/dorado
# Update `dorado_binary:` in pipeline_config.yaml to match.

mamba create -n nanopore $CHANNELS --yes \
    chopper flye minimap2 samtools

# === 2. autocycler — autocycler + 8 assemblers ===
mamba create -n autocycler $CHANNELS --yes \
    autocycler canu flye metamdbg miniasm necat nextdenovo plassembler raven

# === 3. prokka — Bakta annotation ===
# Env name is 'prokka' for back-compat; the actual tool installed is Bakta.
mamba create -n prokka $CHANNELS --yes \
    bakta

# Bakta database (~50 GB; pick `light` to save disk if you don't need all dbs)
mamba run -n prokka bakta_db download \
    --output "$DB_DIR/bakta" --type full
# Update `bakta_db:` in pipeline_config.yaml to:
#   $DB_DIR/bakta/db

# === 4. gtdbtk3 — GTDB-Tk classification ===
mamba create -n gtdbtk3 $CHANNELS --yes \
    gtdbtk
# GTDB-Tk database (~110 GB, R220 release as of this writing).
# The env's `download-db.sh` helper handles it; or set GTDBTK_DATA_PATH manually.
mamba run -n gtdbtk3 download-db.sh "$DB_DIR/gtdbtk"
# Either set GTDBTK_DATA_PATH in the conda env activation, or put the path
# into pipeline_config.yaml under `gtdbtk_data_path:`.

# === 5. mmseqs — minimap2, samtools, metabat2, plotting stack ===
# Everything for s13_map_reads, s15_bin, s17_checkm2 input prep, and s20_report.
mamba create -n mmseqs $CHANNELS --yes \
    minimap2 samtools metabat2 \
    pandas matplotlib numpy

# === 6. checkm2 — genome quality assessment ===
mamba create -n checkm2 $CHANNELS --yes \
    checkm2
# CheckM2 diamond db (~3 GB)
mamba run -n checkm2 checkm2 database \
    --download --path "$DB_DIR/checkm2"
# Set `checkm2_db:` in pipeline_config.yaml to:
#   $DB_DIR/checkm2/CheckM2_database/uniref100.KO.1.dmnd

# === 7. mijamp — methylation analysis (optional) ===
mamba create -n mijamp $CHANNELS --yes \
    pandas meme ont-modkit minimap2 samtools biopython xmltodict
git clone https://code.ornl.gov/alexander-public/mijamp.git "$SOFTWARE_DIR/mijamp"
chmod +x "$SOFTWARE_DIR/mijamp/scripts/preprocess"
# Set `mijamp_dir:` in pipeline_config.yaml to:
#   $SOFTWARE_DIR/mijamp

# === Verify everything is callable ===
for env in nanopore autocycler prokka gtdbtk3 mmseqs checkm2 mijamp; do
    if mamba env list | grep -q "^$env "; then
        echo "[OK]  env $env exists"
    else
        echo "[ERR] env $env not created"
    fi
done
```

Notes / caveats:

- Dorado is a vendor binary; the URL above is for x86_64 Linux. For ARM
  or macOS, grab the right archive from <https://github.com/nanoporetech/dorado>.
- The basecall models referenced in `pipeline_config.yaml` (`dna_r10.4.1_*`)
  must also be downloaded — use `dorado download` after extracting:
  `<dorado_bin> download --model dna_r10.4.1_e8.2_400bps_sup@v5.2.0`.
  The pipeline's `basecall_model` and `mod_models` paths must point to
  wherever you put them.
- `nextdenovo`, `plassembler`, and `metamdbg` are occasionally flaky on
  bioconda — if `mamba create` fails for the `autocycler` env, install the
  difficult tool into its own env or use the upstream install instructions.
- The default `gtdbtk` db is huge. If your `/scratch` quota is tight, use
  the lite db (`download-db.sh --lite`) at the cost of less-precise species
  assignments.
- `mmseqs` is the env name we use, but it's a misnomer — we install
  minimap2 + samtools + metabat2 + plotting libs in it, not mmseqs2 itself.
  Rename in `pipeline_config.yaml` if that bugs you.

## Sample sheet

TSV with these columns (header required, order doesn't matter):

**Required:**

| column              | meaning                                              |
|---------------------|------------------------------------------------------|
| `barcode`           | dorado demux barcode tag, e.g. `barcode01`           |
| `sample`            | output prefix; must be unique                        |
| `genome_length`     | expected genome size in bp (integer; commas allowed) |
| `genus`             | Bakta `--genus`                                      |
| `species`           | Bakta `--species`                                    |
| `strain`            | Bakta `--strain` (also used as `--prefix`)           |

**Optional** (any of these may be absent from the header, or left blank per row):

| column              | meaning                                                            |
|---------------------|--------------------------------------------------------------------|
| `kingdom`           | not used by Bakta; kept for back-compat with old samplesheets       |
| `reference_proteins`| path to a `.gbk` reference; emits `--proteins <path>` when present |
| `gram`              | `pos` / `neg` / `unknown` — maps to Bakta `--gram +/-/?`           |
| `plasmid`           | plasmid name — emits `--plasmid <name>` when present               |

## Running

```bash
# 1. Edit pipeline_config.yaml — confirm dorado path, model, kit, conda
#    env names, prokka_db_dir, email.

# 2. Fill out a samplesheet (copy samplesheet_template.tsv and edit).

# 3. Generate submit scripts.
python3 generate_pipeline.py \
    --pod5-dir /ebio/abt3_scratch/jmarsh/nanopore/<run_id>/pod5 \
    --samplesheet /path/to/samplesheet.tsv \
    --output-dir /ebio/abt3_scratch/jmarsh/nanopore/<run_name> \
    --config pipeline_config.yaml

# 4. Submit. Every stage holds on the previous via -hold_jid.
bash <output-dir>/submit_all.sh

# 5. Watch progress.
tail -f <output-dir>/pipeline.log
```

## Pipeline log

Every stage writes `START / OK / SKIP / WARN / FAIL / DONE` events to a
single file at `<output-dir>/pipeline.log`. Each line has an ISO-8601 UTC
timestamp, a level, the stage name, and a message. Example:

```
2026-05-11T08:14:02Z [START] s05_subsample
2026-05-11T08:14:02Z [SKIP] s05_subsample G0224_i1 (4 subsamples already exist)
2026-05-11T08:14:55Z [OK]   s05_subsample V123_i2 subsampled (genome_size=1800000)
2026-05-11T08:15:30Z [WARN] s05_subsample Mmobile: subsample rc=1
2026-05-11T08:15:30Z [DONE] s05_subsample stage finished with 1 per-sample failures
```

Per-sample failures don't abort the stage — they're logged as `WARN` and
the stage exits 0 so downstream stages can still run on the samples that
succeeded. Stage-level crashes (the script itself dying) get an `[ERROR]`
trap-fired `FAIL` event.

Each stage also greps its SGE stdout log for `error|warn|fail|traceback`
patterns at the end and emits the first 50 hits as `[SCAN]` events, so
unusual lines in tool output bubble up to `pipeline.log` automatically.

## Resume

Every stage checks for existing outputs and skips work that's already
done. Re-running `bash submit_all.sh` after a partial failure picks up
where the previous run stopped.

What each stage looks at to decide whether to skip:

| stage             | skip condition                                                                  |
|-------------------|---------------------------------------------------------------------------------|
| `s01_basecall`    | `basecalls.bam` exists and is non-empty                                         |
| `s02_demux`       | any `*.fastq` / `*.fastq.gz` in `fastq_raw/`                                    |
| `s03_rename`      | `fastq/<sample>.fastq` exists                                                   |
| `s04_filter`      | `fastq_filtered/<sample>.fastq` exists                                          |
| `s05_subsample`   | `autocycler/<sample>/sample_04.fastq` exists (all four subsamples present)      |
| `s06_assemble_*`  | `autocycler/assemblies/<sample>_<asm>_<sub>.fasta` exists per command           |
| `s07_organize`    | `assemblies_fasta/<sample>/<file>.fasta` exists per file                        |
| `s08_compress`    | `autocycler/<sample>/input_assemblies.fasta` exists                             |
| `s09_cluster`     | `autocycler/<sample>/clustering/` exists                                        |
| `s10_trim_resolve`| `cluster_*/5_final.gfa` exists (resolve); trimmed marker (trim)                 |
| `s11_combine`     | `autocycler/<sample>/consensus_assembly.fasta` exists                           |
| `s12_collect`     | `consensus/<sample>.fasta` exists                                               |
| `s13_map_reads`   | `mapping/<sample>/flagstat.txt` + `unmapped.fastq` exist                        |
| `s14_meta_assemble` | `meta_flye/<sample>/assembly.fasta` OR `.skipped_low_reads` sentinel          |
| `s15_bin`         | `bins/<sample>/.binning_done` sentinel                                          |
| `s16_annotate`    | `bakta/<id>/<id>.gff3` exists (per consensus AND per metabat2 bin)              |
| `s17_checkm2`     | `checkm2/output/quality_report.tsv` exists                                      |
| `s18_classify`    | any `gtdbtk/output/gtdbtk.*.summary.tsv` exists (consensuses + bins)            |
| `s19_aggregate`   | `reports/metrics.tsv` always regenerated (cheap; no skip guard)                 |
| `s20_report`      | `reports/report.html` regenerated each run                                      |
| `s02b_demux_bam`  | any `bams/*.bam` exists                                                         |
| `s21_mijamp`      | `mijamp/<sample>/` exists and non-empty                                         |

To **force** re-run for one sample, delete its output for the stage you
want to redo. For example, to retry compress + everything downstream for
`Mmobile`:

```bash
rm -rf <output-dir>/autocycler/Mmobile/input_assemblies.fasta \
       <output-dir>/autocycler/Mmobile/clustering \
       <output-dir>/autocycler/Mmobile/consensus_assembly.fasta \
       <output-dir>/consensus/Mmobile.fasta \
       <output-dir>/prokka/Mmobile
bash <output-dir>/submit_all.sh
```

The `qsub` calls for stages whose samples are all already done will run
through quickly and just log `SKIP` events.

## Tuning the contig filter (compress threshold)

`autocycler compress` rejects an input set where the mean contigs per
assembly exceeds 25. Fragmented inputs from miniasm / raven / canu on
noisy reads routinely break this. `s07_organize` filters every assembly
to drop contigs shorter than `min_contig_length` (default 1000 bp) before
compress sees them.

If you still hit the error, two knobs in `pipeline_config.yaml`:

- `min_contig_length` — raise this (e.g. 5000 or 10000). Removes more
  short fragments. Risk: small plasmids get dropped.
- `max_contigs` — raise above 25. The flag is only emitted when changed
  from 25 (autocycler's default), so older versions of autocycler that
  don't accept the flag still work at the default.

After changing either, regenerate and resubmit. Already-organized
samples will reorganize against the new threshold because the python
filter respects the existing output and only rewrites missing files —
to apply a stricter threshold to already-filtered samples, delete the
`assemblies_fasta/` subdir for those samples first.

## Output layout

```
<output-dir>/
├── pipeline.log                   # unified event log
├── basecalled/basecalls.bam       # s01_basecall
├── fastq_raw/                     # s02_demux
├── fastq/<sample>.fastq           # s03_rename
├── fastq_filtered/<sample>.fastq  # s04_filter
├── autocycler/<sample>/           # s05_subsample, s08_compress, s09_cluster
├── autocycler/assemblies/         # s06_assemble (flat)
├── autocycler/assemblies_fasta/   # s07_organize (per-sample, length-filtered)
├── autocycler/cluster/<sample>/   # s09_cluster (symlink)
├── consensus/<sample>.fasta       # s12_collect
├── bakta/<sample>/                # s13_annotate (Bakta output)
├── mapping/<sample>/              # s14_map_reads — mapped.bam, flagstat.txt, unmapped.fastq
├── meta_flye/<sample>/            # s15_meta_assemble — Flye --meta assembly.fasta
├── bins/<sample>/                 # s16_bin — <sample>_bin.N.fa MetaBAT2 outputs
├── checkm2/output/                # s17_checkm2 — quality_report.tsv
├── gtdbtk/
│   ├── batchfile.tsv              # consensuses + bins
│   └── output/                    # gtdbtk.bac120.summary.tsv, ar53, classify/, etc.
├── reports/
│   ├── metrics.tsv                # s19_aggregate — long-format per-assembly metrics
│   ├── plots/*.png                # s20_report — publication-ready figures
│   └── report.html                # s20_report — final HTML report
├── submit_scripts/                # generated SGE scripts
├── stdout/                        # SGE per-stage stdout
└── submit_all.sh
```

## Stage order and hold dependencies

```
s01_basecall (GPU)
  └─ s02_demux (GPU)
      └─ s03_rename
          └─ s04_filter (chopper)
              └─ s05_subsample
                  └─ s06_assemble_batch_001..NNN  (parallel)
                      └─ s07_organize  (with contig filter)
                          └─ s08_compress
                              └─ s09_cluster
                                  └─ s10_trim_resolve
                                      └─ s11_combine
                                          └─ s12_collect
                                              └─ s13_map_reads (minimap2 → consensus)
                                          # (s02b_demux_bam runs in parallel off s01, feeds s21_mijamp)
                                              └─ s14_meta_assemble (Flye --meta on unmapped)
                                                  └─ s15_bin (MetaBAT2)
                                                      └─ s16_annotate (Bakta on consensuses + bins)
                                                          └─ s17_checkm2
                                                              └─ s18_classify (GTDB-Tk, consensuses + bins)
                                                                  └─ s19_aggregate
                                                                      └─ s20_report (plots + HTML)
                                                                          └─ s21_mijamp (needs s20 done AND s02b_demux_bam done)
```

## Assumptions

- The conda envs named in `pipeline_config.yaml` exist and contain the
  expected tools: `nanopore` → dorado + chopper + flye, `autocycler` →
  autocycler plus every assembler in `assemblers:`, `prokka` → Bakta,
  `gtdbtk3` → GTDB-Tk (with `GTDBTK_DATA_PATH` set), `mmseqs` → minimap2 +
  samtools + metabat2 + pandas/matplotlib, `checkm2` → CheckM2, and (if
  using MIJAMP) `mijamp` → pandas, meme, ont-modkit, minimap2, samtools,
  biopython, xmltodict.

## MIJAMP setup

MIJAMP is optional but enabled when `mijamp_dir` in
`pipeline_config.yaml` points to a valid clone of the MIJAMP repo. Install
once on the cluster:

```bash
git clone https://code.ornl.gov/alexander-public/mijamp.git /path/to/mijamp
mamba create -n mijamp -c bioconda -c conda-forge \
    pandas meme ont-modkit minimap2 samtools biopython xmltodict
```

Then set in `pipeline_config.yaml`:

```yaml
mijamp_dir: /path/to/mijamp
conda_env_mijamp: mijamp
```

Leaving `mijamp_dir` empty makes `s21_mijamp` log a `WARN` and exit
cleanly, so the rest of the pipeline still produces all other outputs.
- Dorado demux output filenames contain the barcode label. The `s03_rename`
  step matches `*<barcode>*.fastq[.gz]`.
- Autocycler subsample produces 4 subsamples named `sample_01.fastq` …
  `sample_04.fastq`.
- SGE job names need a leading letter — stages are named `s01_*` … `s13_*`
  to satisfy this.
