# Nanopore assembly + annotation pipeline

End-to-end SGE pipeline that takes raw POD5 files and a sample sheet, runs
basecalling → demux → QC → Autocycler 8-assembler consensus → Bakta →
read-mapping → unmapped-read meta-assembly (Flye --meta) → MetaBAT2
binning → CheckM2 → GTDB-Tk → structural-heterogeneity genotyping (two or
more versions of a genome in one sample) → small variants between samples and
against reference genomes → minority variants within samples, and produces an
HTML report with publication-ready plots and per-sample summary tables.
Optional consensus polishing (dorado polish) sits between assembly and mapping;
it is off by default because it did not improve this pipeline's consensus
genomes in validation.

## Files

- `generate_pipeline.py` — emits all SGE submit scripts for a given run
- `pipeline_config.yaml` — paths, models, conda envs, SGE resources, filter thresholds
- `samplesheet_template.tsv` — fill in one row per barcode
- `build_type_strain_db.py` — run on an internet-connected machine to fetch
  RefSeq type-strain proteomes for `--proteins` (see "Type-strain proteomes")
- `sv_genotype.py` — structural-heterogeneity engine run by `s18b_sv` (copied
  into each run's `submit_scripts/`; also usable on its own, see
  "Structural heterogeneity")
- `polish_tools.py` — read subsampling, read groups, change tables and install
  steps for `s12b_polish` (see "Consensus polishing")
- `snv_tools.py` — comparison planning, calling, annotation and run tables for
  `s18c_snv` (see "Small variants between samples and references")
- `minor_variants.py` — minority-variant engine run by `s18d_minor` (see
  "Minority variants within samples")
- `snv_validation/` — scripts, pre-specified decision rules and result tables of
  the validations that set the defaults of these three stages
- `basecall_split_validation/` — scripts, the pre-specified rule and result tables
  of the validation of split basecalling (see "Split basecalling across GPU nodes")

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

# === 8. sv — structural heterogeneity (s18b_sv) ===
mamba create -n sv $CHANNELS --yes python=3.12 sniffles mappy pysam samtools

# === 9. snv — polishing, small variants and minority variants (s12b, s18c, s18d) ===
mamba create -n snv $CHANNELS --yes python=3.12 pysam numpy scipy pandas mappy \
    minimap2 samtools bcftools htslib vcfdist

# === 10. clair3 — Clair3 >= 2.0 (PyTorch) for s18c ===
mamba create -n clair3 $CHANNELS --yes "clair3>=2.0.2"
# R10.4.1 sup v5.2.0 model (Rerio model converted to PyTorch by the Clair3 authors):
M=$SOFTWARE_DIR/clair3_models/r1041_e82_400bps_sup_v520; mkdir -p $M
URL=https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/r1041_e82_400bps_sup_v520/
for f in pileup.pt full_alignment.pt; do curl -fsSL -o $M/$f $URL$f; done
# Set `snv_clair3_model:` in pipeline_config.yaml to $M. The dorado polishing models
# are looked up in `polish_models_dir:` (default: the folder holding basecall_model);
# fetch them with `<dorado_bin> download --model dna_r10.4.1_e8.2_400bps_polish_bacterial_methylation_v5.0.0`.

# === 11. phava — phase variation (s18e_phava) ===
# (conda-forge before bioconda, as this env was built and tested; PhaVa's own test: `phava test`)
mamba create -n phava -c conda-forge -c bioconda --yes "phava=0.2.3"

# === Verify everything is callable ===
for env in nanopore autocycler prokka gtdbtk3 mmseqs checkm2 mijamp sv snv clair3 phava; do
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
| `reference_proteins`| your own `--proteins` reference — see below                        |
| `gram`              | `pos` / `neg` / `unknown` — maps to Bakta `--gram +/-/?`           |
| `plasmid`           | plasmid name — emits `--plasmid <name>` when present               |

### Supplying your own `--proteins` reference

Put an absolute path in `reference_proteins` and it is passed straight to
whichever annotator runs for that sample:

- **Amino-acid FASTA (`.faa`) and GenBank (`.gbk`, `.gbff`, `.gb`) both
  work**, for both Bakta and Prokka — the same column covers bacterial and
  archaeal samples, so you do not need to know in advance which tool will be
  chosen.
- **Plain FASTA headers are fine; the pipeline reformats them if needed.**
  Bakta requires headers shaped `>SeqID gene~~~product~~~dbxrefs` (a space
  after the ID, then three `~~~`-separated fields) and silently discards
  every record that isn't (`skipped incorrectly formatted user proteins: N,
  imported: 0`, followed by a diamond crash). A normal RefSeq header such as
  `>WP_000123456.1 chorismate mutase [Escherichia coli]` therefore has to be
  rewritten. The annotation stage detects this and writes a converted copy to
  `<output>/bakta/proteins_prepared/`, which is cached and reused. Anything
  already in Bakta's format is used as-is, and Prokka gets the file as-is
  because it tolerates plain FASTA.
- **GenBank is converted for Bakta and passed straight to Prokka.** Bakta 1.12
  reads `--proteins` as FASTA first and exits if that fails, and Biopython
  ≥ 1.87 fails on any GenBank file, so a GenBank reference handed directly to
  Bakta dies with `User proteins file Fasta format not valid!` (the stage then
  re-annotates without `--proteins` and logs a `WARN`). For Bakta the
  annotation stage therefore extracts the CDS translations into
  `<output>/bakta/proteins_prepared/<name>.<hash>.bakta.faa`, following the
  same rules as Bakta's own GenBank reader: `/pseudo` and `/pseudogene` CDSs
  are skipped, `/locus_tag`, `/product` and `/translation` are required, and
  `/db_xref` plus `/EC_number` become dbxrefs. The copy is cached per source
  file (path, size and modification time). Prokka converts GenBank itself.
  EMBL files are still passed through unchanged.
- It **overrides** the GTDB-derived RefSeq type-strain proteome. Leave the
  cell blank to let the pipeline pick the type strain automatically.
- The path is on the **cluster**, and is checked at run time. If the file is
  missing the stage logs a warning and annotates without `--proteins` rather
  than failing.

This column is keyed by sample, so it only covers the Autocycler
consensuses. MetaBAT2 bins do not exist until the pipeline has run, so they
always take the automatic type-strain lookup.

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

At the end, each stage also greps its SGE stdout log for definitive error
markers (`ERROR:`, `Error:`, `FATAL:`, `FAIL:`, `panic:`, segfaults, `Killed`,
minus known benign assembler chatter) and emits the first 20 hits as `[SCAN]`
events, so unusual lines in tool output bubble up to `pipeline.log`
automatically. SGE appends every run of a stage to the same
`stdout/<stage>.log`, so only the lines written by the current run are
scanned; the line number in each `[SCAN]` event refers to the whole file.

`bash <output-dir>/pipeline_status.sh` reports the latest attempt of each
stage (its events from its last `[START]` onwards): stages started but not
yet `DONE` (still running, failed, or killed by SGE), stages whose `DONE`
reports per-sample failures, stage-level `FAIL`s, and per-level event
counts. Events from earlier attempts stay in `pipeline.log` as history but
are not reported again once the stage has been re-run. At the end of each
run `s20_report` writes the same summary to `<output-dir>/pipeline_summary.txt`.
Logs written by earlier versions of the pipeline contain `pipeline summary`
blocks (copies of past `START`/`DONE`/`FAIL` events that `s20_report`
appended to `pipeline.log`); `pipeline_status.sh` ignores them.

## Resume

Every stage checks for existing outputs and skips work that's already
done. Re-running `bash submit_all.sh` after a partial failure picks up
where the previous run stopped.

What each stage looks at to decide whether to skip:

| stage             | skip condition                                                                  |
|-------------------|---------------------------------------------------------------------------------|
| `s01_basecall`    | `basecalls.bam` exists and is non-empty                                         |
| `s01_basecall` (split, task i) | `basecalls.bam` exists, or `basecalled/chunks/chunk_<i>.bam` with a `.done` record for the same chunk files and basecalling settings |
| `s01b_merge_bam`  | `basecalls.bam` exists and is non-empty (written only after every chunk has been checked) |
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
| `s16_classify`    | any `gtdbtk/output/gtdbtk.*.summary.tsv` exists (consensuses + bins)            |
| `s17_checkm2`     | `checkm2/output/quality_report.tsv` exists                                      |
| `s18_annotate`    | `bakta/<id>/<id>.gff3` (Bakta) or `.gff` (Prokka), per genome                   |
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
├── pipeline_summary.txt           # s20_report — latest attempt of each stage
├── basecalled/basecalls.bam       # s01_basecall (split basecalling: s01b_merge_bam)
├── basecalled/chunks/             # split basecalling: chunk_<i>.done, summary.tsv, merged_header.sam
│                                  #   (chunk_<i>.bam until s01b_merge_bam has merged them)
├── fastq_raw/                     # s02_demux
├── bams/                          # s02b_demux_bam — per-barcode BAMs (keep MM/ML tags)
├── fastq/<sample>.fastq           # s03_rename
├── fastq_filtered/<sample>.fastq  # s04_filter
├── autocycler/<sample>/           # s05_subsample, s08_compress, s09_cluster
├── autocycler/assemblies/         # s06_assemble (flat)
├── autocycler/assemblies_fasta/   # s07_organize (per-sample, length-filtered)
├── autocycler/cluster/<sample>/   # s09_cluster (symlink)
├── consensus/<sample>.fasta       # s12_collect (replaced by the polished genome if s12b installs it)
├── polish/<sample>/               # s12b_polish — draft.fasta, polished.fasta, changes.tsv, variants.vcf, info.json
├── mapping/<sample>/              # s13_map_reads — mapped.bam, flagstat.txt, unmapped.fastq
├── meta_flye/<sample>/            # s14_meta_assemble — Flye --meta assembly.fasta
├── bins/<sample>/                 # s15_bin — <sample>_bin.N.fa MetaBAT2 outputs
├── gtdbtk/                        # s16_classify
│   ├── batchfile.tsv              # consensuses + bins
│   └── output/                    # gtdbtk.bac120.summary.tsv, ar53, classify/, etc.
├── checkm2/output/                # s17_checkm2 — quality_report.tsv
├── bakta/                         # s18_annotate
│   ├── <genome_id>/               # Bakta (Bacteria) or Prokka (Archaea) output
│   ├── annotation_plan.tsv        # annotator and --proteins chosen per genome
│   └── proteins_prepared/         # --proteins references converted for Bakta (cached)
├── reports/
│   ├── metrics.tsv                # s19_aggregate — long-format per-assembly metrics
│   ├── plots/*.png                # s20_report — publication-ready figures
│   └── report.html                # s20_report — final HTML report
├── sv/                            # s18b_sv — structural heterogeneity
│   ├── events.tsv                 # one row per sample x locus: allele fractions, status, segment, genes
│   ├── sample_summary.tsv         # per-sample counts by status + reported events
│   ├── null_loci.tsv              # fake events used for each sample's detection floor
│   ├── summary.json               # everything above, machine-readable
│   ├── sv_genotype.log
│   └── <sample>/                  # sniffles.vcf, alleles.fasta, reads/<locus>.<allele>.txt
├── snv/                           # s18c_snv — small variants between samples and references
│   ├── comparisons.tsv            # one row per comparison (self / cross / reference)
│   ├── calls.tsv                  # every call: read support, repeat/homopolymer flags, gene effect
│   ├── plan.tsv, similarity.tsv   # chosen comparisons and all measured genome identities
│   └── <comparison>/              # calls.vcf, calls.tsv, assembly_diffs.tsv, aln.bam, ref.fasta
├── minor/                         # s18d_minor — minority variants within samples
│   ├── summary.tsv                # one row per sample: positions assessed, depth, detection floor, counts
│   └── <sample>/                  # sites.tsv, groups.tsv, summary.tsv, null_model.tsv
├── phava/                         # s18e_phava — phase variation (invertons)
│   ├── summary.tsv                # one row per sample: candidates, spanned, called by class
│   ├── invertons_called.tsv       # called invertons of all samples, linked across samples
│   ├── invertons_all.tsv          # every candidate of every sample
│   └── <sample>/                  # invertons.tsv, summary.tsv, phava_ratio.tsv, phava_IRs.tsv, phava_genes.gff
├── mijamp/<sample>/               # s21_mijamp — MIJAMP preprocessing (mapped.bam, data.bed)
├── submit_scripts/                # generated SGE scripts
├── stdout/                        # SGE per-stage stdout
├── submit_all.sh
└── pipeline_status.sh             # status of the latest attempt of each stage
```

## Stage order and hold dependencies

```
s01_basecall (GPU; basecall_split_jobs > 1: an array, one GPU job per chunk of pod5 files)
  └─ [s01b_merge_bam]  (split only: checks the chunks and concatenates them; s02 and s02b then wait on it)
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
                                              └─ [s12b_polish]  (optional, polish_enabled; s13 then waits on it)
                                              └─ s13_map_reads (minimap2 → consensus)
                                          # (s02b_demux_bam runs in parallel off s01 or s01b, feeds s21_mijamp)
                                              └─ s14_meta_assemble (Flye --meta on unmapped)
                                                  └─ s15_bin (MetaBAT2)
                                                      ├─ s16_classify (GTDB-Tk, consensuses + bins)
                                                      │      └─ s18_annotate (Bakta / Prokka)
                                                      │             ├─ s18b_sv (needs s13 mapping + s18 annotation)
                                                      │             ├─ s18c_snv (optional; needs s13 + s18)
                                                      │             ├─ s18d_minor (optional; needs s13 + s18)
                                                      │             └─ s18e_phava (optional; needs s13 + s18)
                                                      └─ s17_checkm2 (consensuses + bins)
                                                                  └─ s19_aggregate (waits on s17 + s18 + s18b [+ s18c + s18d + s18e])
                                                                      └─ s20_report (plots + HTML)
                                                                          └─ s21_mijamp (needs s20 done AND s02b_demux_bam done)

GTDB-Tk and CheckM2 both take the raw assemblies and neither needs the other,
so they run in parallel. Annotation genuinely depends on classification —
GTDB-Tk's call decides both which annotator runs and which reference proteome
is passed to `--proteins`.

## Split basecalling across GPU nodes (s01_basecall + s01b_merge_bam)

dorado runs on one node, so a single `s01_basecall` job works through the whole run on whatever card
it gets. With `basecall_split_jobs: N` (N > 1), `generate_pipeline.py` divides the pod5 files into N
chunks of similar size, `s01_basecall` becomes an SGE array with one GPU job per chunk, and a new
CPU stage, `s01b_merge_bam`, concatenates the chunk BAMs into `basecalled/basecalls.bam`. Demux and
everything after it are unchanged. The `gpu` complex is 1 per host here, so each job runs on its own
node; as many run at once as there are free GPU nodes, and the rest queue.

1. **Plan, at generate time.** The `*.pod5` files directly inside `--pod5-dir` (the files dorado
   reads; it does not recurse) go largest first to the chunk with the fewest bytes so far, so chunk
   sizes differ by at most one file. The plan is written to `submit_scripts/s01_basecall_chunks.tsv`
   (chunk, file, bytes) and printed; for 20261012_cm_conjugation:
   `s01_basecall: 32 GPU jobs over 83 pod5 files (455.8 GB); chunks of 11.7-17.5 GB`.
2. **Chunk jobs (`s01_basecall`, array task i).** Task i links chunk i's files into
   `basecalled/chunks/pod5_<i>/` and stops if one is missing or has changed size. It runs the
   single job's dorado command (same model, modification models and kit) on that folder into
   `chunk_<i>.bam.partial`, which becomes `chunk_<i>.bam` only when dorado exits 0. Then
   `chunk_<i>.done` records the files, bytes, host, GPUs, run time and a key of the chunk's files and
   the basecalling settings. Each task logs to `pipeline.log` as `s01_basecall.<i>`, so
   `pipeline_status.sh` lists the chunks still running; its SGE log is `stdout/s01_basecall.<i>.log`.
3. **Merge (`s01b_merge_bam`, CPU, holds on `s01_basecall`).** It checks that the pod5 folder holds
   exactly the planned files and that every chunk finished under the current plan and settings and
   passes `samtools quickcheck -u`. It builds the header from chunk 1's `@HD` and `@PG`, the `@RG`
   lines of all chunks (an id must carry the same line in every chunk) and `@CO` lines that record the
   split. It concatenates the chunks in order with `samtools cat`, checks that the record count
   equals the chunks' total, and only then moves the result to `basecalls.bam`. Per-chunk reads, times
   and hosts go to `basecalled/chunks/summary.tsv`. The chunk BAMs duplicate `basecalls.bam` and are
   then deleted, unless `basecall_split_keep_chunks: true`.

| key | shipped | when absent | meaning |
|---|---|---|---|
| `basecall_split_jobs` | `32` | `1` | number of GPU jobs; `1` = one job over the whole folder and no `s01b` |
| `basecall_split_gpumem` | `"12G"` | `sge_gpumem` | `gpumem` request of each job; 12G admits the 12 GB cards as well as the 24 and 80 GB ones |
| `basecall_split_h_rt` | `"24:0:0"` | `sge_h_rt_gpu` | run-time limit of each job; up to 24 h admits standard.q as well as long.q |
| `basecall_split_slots` | `4` | `sge_pe_parallel` | slots of each job; dorado keeps ~4.5 cores busy, and the 12 GB nodes rarely have more free |
| `basecall_split_h_vmem` | `"80G"` | `sge_h_vmem_gpu` | per slot; must stay large (below) |
| `basecall_split_keep_chunks` | `false` | `false` | keep `basecalled/chunks/chunk_*.bam` after the merge |
| `sge_h_vmem_merge_bam` | `"4G"` | `"4G"` | per slot of `s01b_merge_bam` (4 slots) |

### GPU nodes and run time

| nodes | queue | cards per node | pod5 per hour (sup v5.2 + 4mC_5mC + 6mA) | one 14.2 GB chunk |
|---|---|---|---|---|
| node514–519 | long.q | 2 × 24 GB | 20.5 GB/h (run 20261012 on 2 × RTX A5000: 455.8 GB in 22.2 h) | 0.7 h |
| node520–523, node701–705 | long.q | 3–4 × 11–12 GB | 6.8 GB/h (run 20261020 on node703, 3 × RTX 2080 Ti: 77.3 GB in 11.4 h) | 2.1 h |
| node601–609 | standard.q (≤ 24 h) | 4 × TITAN Xp 12 GB | 2.9 GB/h (measured in the validation) | 4.9 h |
| node527, node528 | standard.q (≤ 24 h) | 1–2 × 80 GB | not measured | |

- For this run's 455.8 GB, `basecall_split_jobs: 32` makes chunks of 11.7–17.5 GB (2–3 pod5
  files each). With enough free GPU nodes the basecalling ends when the slowest chunk does: the largest
  chunk takes about 0.9 h on a 24 GB node, 2.6 h on an 11–12 GB long.q node and 6.0 h on a
  TITAN Xp node. The production run's single job took 22.2 h on a 24 GB node; on an RTX 2080 Ti node
  one job would take about 67 h.
- On the validation data (6.47 GB, one TITAN Xp node class), one job took 2 h 13 min
  (2.9 GB/h), three jobs 52 min (chunks of 1.7–2.5 GB), and the merge 3 s.
- More chunks than free nodes is fine: the remaining jobs start as nodes free up, and smaller chunks
  keep a chunk on a slow card from holding up the run. The generator warns when the largest chunk
  could take longer than `basecall_split_h_rt` at the slowest measured speed (2.8 GB/h).
- `h_vmem` is enforced as each job's address-space limit, and dorado's CUDA start-up needs far more
  address space than the at most ~26 GB it then uses. At 4 × 10G (40 GB), dorado stopped within
  seconds with "CUDA out of memory" on a TITAN Xp node and on an 11 GB-card node; at 4 × 80G it ran.
  The per-host `h_vmem` pools hold terabytes, so the large request does not delay scheduling.
- To keep fewer GPU nodes busy at once, cap the running tasks with `qalter -tc 8 <job id of s01_basecall>`.

### Resume and failures

| situation | what happens |
|---|---|
| a job is killed (h_rt, node failure, qdel) | it leaves only `chunk_<i>.bam.partial`; the next `submit_all.sh` basecalls that chunk again from the start |
| a pod5 file is missing or has changed size | its job stops before dorado starts (WARN and FAIL in `pipeline.log`) |
| a chunk did not finish | `s01b_merge_bam` names it and stops; `basecalls.bam` is not written |
| a chunk BAM is truncated or corrupt | the merge moves it to `chunk_<i>.bam.bad`, removes its `.done` and stops; the next `submit_all.sh` basecalls it again |
| the pod5 folder has files the plan lacks (for example, the sequencing run was still writing them) | the merge stops, since their reads would be missing. Re-run `generate_pipeline.py` with the same arguments: it keeps the finished chunks and plans only the new files into extra chunks, so the next `submit_all.sh` basecalls just those |
| `submit_all.sh` after a finished run | every job and the merge skip, since `basecalls.bam` exists |

To start over, delete `basecalled/` and `submit_scripts/s01_basecall_chunks.tsv`. While the plan's
files are unchanged, re-generating keeps the plan.

### Validation

Rule R5 (`basecall_split_validation/RULES.md`) was fixed before the runs, with two dated amendments:
A1 sets the shipped resources, A2 raises `h_vmem` after the first submission failed as described
above. It was tested on the four smallest pod5 files of 20261012_cm_conjugation (6.47 GB). Arm
S3 (3 chunks + merge), arm S1 (the single job) and, because some records differed, a second
single-job arm S1′ ran on the same node class (TITAN Xp).

| check | result |
|---|---|
| reads | the same 36960 reads in S3 and S1, each once (0 only in either) |
| records (FLAG, SEQ, QUAL, every tag incl. MM/ML, RG, barcode) | 36916 of 36960 identical; the other 44 differ only in the modification tags: ML in 44 reads (in 42 of them 1–7 per-base probabilities changed, by at most 7 of 255) and MM in 2. The single job run again (S1′) differs from S1 in the same 44 reads and from S3 in none: all 36960 S3 records equal S1′'s. These differences are run-to-run variation of the single job, not the split; R5 criterion 2 (at most 63.9) is met |
| merged file | `samtools quickcheck -u` passes; 36960 records = the 3 chunks' total = S1; all 6 read groups used are declared (the same 25 `@RG` lines as S1); `@PG` ids unique; 4 `@CO` lines |
| demux on the merged BAM | `s02_demux`: 6 FASTQ files with 36960 reads, and `s02b_demux_bam`: 6 BAMs with 36960 reads; the same read ids per file as from S1 |
| failure handling | a missing chunk, a truncated chunk and a pod5 file outside the plan each stopped the merge (FAIL, no `basecalls.bam`); a changed file size stopped its job before dorado; after the truncated chunk, re-submitting re-basecalled only that chunk and restored S1's 36960 reads; after a fifth file was added, re-generating kept chunks 01–03 and re-submitting basecalled only chunk 04, giving S1's 36960 reads plus the fifth file's 13254 (the same parent reads as production); re-submitting the finished run skipped every job and the merge |
| against the production single-job BAM (other hardware) | R5 asks for the same parent reads, and this is not met: 1 parent read is only in S3 (35 bases, mean Q 3.98) and 1 only in production (177 bases, Q 7.45). S1 differs from production in exactly the same reads, parent reads and fields, so the GPU hardware causes this, not the split. Across hardware none of the 36833 shared reads has an identical record (SEQ differs in 35020, QUAL in 36830) |

Scripts, the rule and the result tables are in `basecall_split_validation/`.

## Annotation: Bakta, Prokka, and --proteins

After GTDB-Tk classifies every assembly (consensuses and MetaBAT2 bins alike),
`s18_annotate` decides per genome:

| GTDB-Tk domain | annotator | notes                                    |
|----------------|-----------|------------------------------------------|
| Bacteria       | Bakta     | plus `--gram` / `--plasmid` if in sheet  |
| Archaea        | Prokka    | `--kingdom Archaea`                      |
| unclassified   | Bakta     | safer general-purpose default            |

Both tools live in the same conda env (`conda_env_bakta`, default `prokka`).

Every genome gets `--proteins <RefSeq type-strain proteome>` for the species
GTDB-Tk assigned. Resolution cascades:

1. an explicit `reference_proteins` in the samplesheet — your own `.faa` or
   `.gbk`, applied whether the genome is routed to Bakta or Prokka (wins)
2. exact species match in the local type-strain index
3. species match after stripping GTDB's alphabetic genus/species suffixes
   (`Methanobrevibacter_A smithii` → `Methanobrevibacter smithii`)
4. any type strain in the same genus
5. nothing — the genome is annotated without `--proteins`, and the species is
   appended to `bakta/missing_species.txt`

Organism metadata (`--genus` / `--species` / `--strain`) comes from the
samplesheet for isolate consensuses and from GTDB-Tk for bins.

### Type-strain proteomes (offline, two-phase)

The cluster has no internet, so the proteomes have to be fetched elsewhere
and copied in. `build_type_strain_db.py` does that.

**Phase 1 — run the pipeline once.** GTDB-Tk classifies everything;
`s18_annotate` annotates without `--proteins` and writes
`<output>/bakta/missing_species.txt`.

**Phase 2 — build the DB on a machine with internet:**

```bash
python3 build_type_strain_db.py \
    --from-gtdbtk <output>/gtdbtk/output/gtdbtk.bac120.summary.tsv \
    --from-gtdbtk <output>/gtdbtk/output/gtdbtk.ar53.summary.tsv \
    --out-dir ./type_strain_db

# or, driven by what the pipeline said it was missing:
python3 build_type_strain_db.py \
    --species-list missing_species.txt \
    --out-dir ./type_strain_db
```

It parses NCBI's `assembly_summary_refseq.txt`, keeps only assemblies flagged
as type material (preferring `assembly from type material` > synonym >
pathotype > neotype, then Complete Genome > Chromosome > Scaffold > Contig),
and downloads each `*_protein.faa.gz`. Repeat runs are additive, so you can
keep topping the collection up as new species show up.

Then copy it to the cluster:

```bash
rsync -av ./type_strain_db/ <cluster>:/ebio/abt3_scratch/$USER/type_strain_db/
```

**You usually don't need to touch the config.** If `type_strain_db_dir` is
empty, the annotation stage searches these locations at runtime and adopts
the first one containing a valid `index.tsv`:

1. `<output-dir>/type_strain_db` (always tried first)
2. everything in `type_strain_db_search_paths`, which defaults to
   `/ebio/abt3_scratch/$USER/type_strain_db`,
   `/ebio/abt3_projects2/Gut_nanopore/data/type_strain_db`, and
   `~/type_strain_db`

An explicit `type_strain_db_dir` still wins when it's usable; if it points
somewhere without an `index.tsv` the stage logs a warning and falls back to
the search list rather than silently annotating without `--proteins`. The
log records which DB was used and how it was found:

```
[OK] s18_annotate type_strain_db: /ebio/abt3_scratch/jmarsh/type_strain_db (auto-discovered existing DB)
```

Then delete `<output>/bakta/` and re-run; annotation picks up the proteomes.
Everything upstream (assembly, classification) skips.

Re-running `build_type_strain_db.py` against an existing directory is
additive and safe: entries already present are reused, entries whose `.faa`
has gone missing are re-fetched, and only genuinely new species are
downloaded.
```

## Structural heterogeneity (s18b_sv)

`s18b_sv` asks, for every sample: is more than one version of this genome
present, and in what proportions? Typical cases are a strain mixed with its
deletion mutant, a subpopulation carrying an inversion, a tandem amplification,
an element that excises and circularises, or an IS element that has moved in
part of the culture. It runs `sv_genotype.py` (copied into
`submit_scripts/_sv_genotype.py`) on each sample's Autocycler consensus and its
s13 read mapping, after s18 so that genes in variable segments can be listed.

### Method

1. **Candidates.** (a) Differences between the genome versions listed in
   `sv_variant_genomes` (single-sequence FASTAs; the first is the base). They
   are found by collinear whole-genome alignment, after rotating and
   reverse-complementing each version to the base. (b) With `sv_discovery`,
   Sniffles2 in mosaic mode (`--mosaic-af-min sv_min_fraction`) on each
   sample's reads. Sniffles2 only proposes candidates; it does not decide
   anything. Calls that span the origin of a circular contig (split reads at
   the contig ends) are skipped. Every locus is placed in every sample's
   consensus through unique 2-kb flank probes, so it is genotyped wherever it
   can be placed, whichever sample it came from. Coordinates are reported on
   each sample's own consensus.
2. **Explicit alleles.** Each locus is a set of full allele sequences: a shared
   left flank, one version of the variable segment per allele, and a shared
   right flank. When a discovered deletion or duplication has its breakpoints
   in a direct repeat pair (IS elements, rRNA operons), the allele that keeps
   one repeat copy at the junction is added next to the allele exactly as
   called. Split alignments in a repeat can overlap in read coordinates, so
   callers often report breakpoints at the outer repeat edges. Reads decide
   between the two, and their counts are summed.
3. **Spanning reads.** Reads at the locus are re-aligned (mappy) to every
   allele. A read counts for an allele only if one alignment covers that
   allele's whole window with ≥ 85 % identity and no internal indel as large as
   the allele difference. Windows are the same length for every allele (the
   shorter allele is padded with flank), so a genome copy carrying either
   version has the same chance of producing a counted read. Both window ends
   carry ≥ `sv_anchor` bp of sequence that is unique in the genome (and within
   each allele), so repeats and tandem copies cannot place a read on the wrong
   allele. A read must also stay collinear with its allele beyond the window:
   a read that leaves the allele inside its flanks (≥ 500 unaligned bp where
   the allele continues) is not counted for any allele. This excludes chimeric
   reads, reads from plasmids carrying homology arms, and reads with another
   difference nearby. Reads are recruited only if their alignments cover both
   anchors, and at most `sv_max_reads` are used (a deterministic subsample).
4. **Junctions, for events longer than the reads.** If fewer than
   `sv_min_spanning` reads span the full windows, or the windows are longer
   than `sv_max_span_window` (20 kb; large deletions, amplifications,
   excisable elements), `sv_junction_mode: auto` makes the breakpoint
   junctions the primary test. Spanning-read counts are still reported. These
   are: left flank|segment start, segment end|right flank, segment
   end|segment start (a tandem copy or a circular form), and left flank|right
   flank (segment absent). For inversions they are the reference and inverted
   junctions at both ends. All junction alleles share one window geometry that
   reaches unique sequence on both sides of any repeat next to the junction.
   Each read region that spans a junction window supports one junction allele.
   Outputs: the fraction of genomes lacking the segment, and end-to-start
   junctions per genome. Junctions alone cannot tell a tandem copy from a
   circle; the coverage ratio of the segment is reported alongside.
5. **Fractions and detection floor.** Allele fractions are shares of genome
   copies, with Wilson 95 % intervals. For each sample and event class (type ×
   size bin), `sv_null_loci` fake events of the same kind and size are
   genotyped at random unique positions of that sample's consensus. The highest
   fake-allele fraction there, but never less than `sv_min_fraction`, is the
   sample's detection floor. An allele is called present when it has ≥
   `sv_min_alt_reads` reads and the lower bound of its interval exceeds the
   floor.

### Outputs and statuses

`sv/events.tsv` has one row per sample × locus: the locus and its source
(`known` / `discovered`), its position in that sample's consensus, which allele
the consensus carries, reads per allele with fractions and intervals, the
floor, the status, flags, genes overlapping the variable segment
(`genes_in_consensus_segment`, from the s18 GFF; a segment that runs across the
origin of a circular contig counts the genes on both sides), the junction-mode
columns, and two columns that place the event among the genes:

| column | content |
|---|---|
| `consensus_segment` | the event's variable sequence in this sample's consensus as `contig:start-end` (1-based, inclusive; start > end runs across the origin), or `contig:p^q` when the consensus lacks that sequence (an insertion carried by part of the population, or a consensus that is the deleted version): the junction between consensus bases p and q where the sequence joins |
| `breakpoint_genes` | for a segment `start contig:S <context> \| end contig:E <context>`, for a junction `junction contig:p^q <context>`. Context: `in locus_tag(gene):product` for each gene containing the site, with ` [cut]` when the gene extends beyond the segment at that end (the breakpoint splits it); otherwise `intergenic (left <gene>, n bp; right <gene>, m bp)`: the nearest gene on each side and the number of bases between it and the site (wrapping around circular contigs) |

Genes are the GFF features other than `gene`, `region` and `source` lines (CDS, tRNA, rRNA, …), the
same features in all three columns. In the report, each reported event shows its segment, the complete
list of genes in the segment (folded when longer than three) and its breakpoint genes.

Statuses:

| status | meaning |
|---|---|
| `mixture` | two or more alleles present above the floor |
| `single:<allele>` | one version only (known loci) |
| `consensus_discrepancy` | discovered locus: reads support only the non-consensus allele (the consensus is probably wrong here, e.g. a collapsed tandem copy) |
| `not_supported` | discovered candidate not confirmed by re-genotyping |
| `insufficient_spanning_reads` / `insufficient_junction_reads` | fewer than `sv_min_spanning` counted reads |
| `repeat_limited` | no unique anchor within reach (e.g. the segment sits inside a repeat longer than the reads) |

The `method` column says whether the status comes from spanning reads or
junctions. With `sv_junction_mode: always` both are computed for every locus.
Flags: `anchor_extended_past_repeat`, `spans_contig_origin`,
`strand_bias:<allele>`, `primary_elsewhere:N` (reads whose primary alignment
lies elsewhere), `divergent_reads:N` (≥ 5 % of reads excluded by the
collinearity check), `consensus_matches_no_allele`. Read names per allele are
in `sv/<sample>/reads/` for inspection in IGV or Geneious.

### Two versions of one strain (e.g. a strain and its deletion mutant)

Put both genome versions in `sv_variant_genomes`. The base is listed first; it
does not need to match the consensus in rotation or strand.

```yaml
sv_variant_genomes:
  - /path/to/T002_i2.fasta
  - /path/to/T002_i2_with_deletion.fasta
```

Each sample then gets a row for the difference, here a 6,527-bp deletion,
with the fraction of genome copies carrying each version. The result does not
depend on which version Autocycler assembled as the consensus: in the
validation below the two choices gave estimates within 0.3 percentage points.

If the mutant was made by allelic exchange, reads from the knock-out
plasmid (or a single-crossover integrant) contain the deletion junction
followed by homology arm and then vector. The collinearity check removes such
reads when the vector starts inside the allele flanks, which reach ≥ 8 kb from
the junction. For extra safety, set `sv_anchor` longer than the homology arms
(e.g. arm length + 500). The `junction` columns use short windows; if they
report a deletion that the spanning columns do not, suspect construct-borne
junctions before a chromosomal subpopulation.

### Running it by hand

```bash
mamba activate sv
printf 'S1\t<consensus.fasta>\t<mapped.bam>\t<annotation.gff>\t<annotation.fna>\n' > samples.tsv
python sv_genotype.py --samples samples.tsv --outdir sv_out --discovery \
    --variant-genomes versionA.fasta versionB.fasta --threads 8
```

The GFF and annotation-FASTA columns are optional; the FASTA maps annotation
contig names to consensus contigs by sequence.

### Validation

Simulated mixtures of *M. smithii* T002 genome versions (Badread 0.4.2, `nanopore2023` error and q-score models,
read identity mean 98 %, sd 2 %, max 99.8 % — lower than this run's reads, median 99.7 % — and mean length
11 kb, this run: 11.5 kb), mapped with minimap2
`-ax map-ont` as in s13. Events: the 6,527-bp deletion of `T002_i2_with_deletion` (genotyped on either version
as the consensus), a synthetic 3-kb inversion and a synthetic 5-kb tandem duplication. Fractions of the
alternative version: 0, 1, 2, 5, 10, 25, 50, 75, 90, 95, 98 and 100 %; read depths 100×, 400× (all events) and
1,000× (deletion). Each condition was run twice: with both genome versions supplied (known locus) and with
discovery only. Truth is the realised share of reads from each version in the mixture.

| event | depth | spanning reads (median) | mean abs. error, 5–95 % (pp) | bias (pp) | 95 % CI contains truth | lowest fraction called | false mixture calls in pure samples | discovered and reported, 5–95 % |
|---|---|---|---|---|---|---|---|---|
| DEL (consensus A) | 100× | 40 | 4.2 | +1.2 | 10/10 | 10% | 0 | 6/7 |
| DEL (consensus A) | 400× | 175 | 0.8 | -0.1 | 10/10 | 5% | 0 | 7/7 |
| DEL (consensus A) | 1000× | 432 | 1.1 | -0.7 | 10/10 | 2% | 0 | 7/7 |
| DEL (consensus B) | 100× | 40 | 4.2 | +1.2 | 10/10 | 10% | 0 | 6/7 |
| DEL (consensus B) | 400× | 175 | 0.9 | -0.2 | 10/10 | 5% | 0 | 7/7 |
| DEL (consensus B) | 1000× | 431 | 1.1 | -0.7 | 10/10 | 2% | 0 | 7/7 |
| DUP | 100× | 31 | 3.7 | +3.4 | 10/10 | 25% | 0 | 5/7 |
| DUP | 400× | 123 | 1.8 | -1.0 | 10/10 | 5% | 0 | 7/7 |
| INV | 100× | 65 | 1.8 | -1.2 | 10/10 | 10% | 0 | 6/7 |
| INV | 400× | 251 | 0.8 | -0.3 | 10/10 | 2% | 0 | 7/7 |

Over all 100 mixtures, the spanning-read 95 % interval contained the true fraction in 100 and the junction interval in 100. No pure sample (20 conditions) was called a mixture, 0 of 21,010 spanning reads were assigned to the wrong version (read names carry the true source), and discovery reported 0 false events. The lowest detectable fraction is set by the number of spanning reads, because an allele needs ≥ 5 reads with its lower bound above the 1 % floor: at 100× from 10 % (deletion, inversion) or 25 % (duplication), at 400× from 2–5 %, and at 1,000× from 2 %. In this run's *M. smithii* samples, about 1,100–1,700 reads span the 6.5-kb deletion window; with that many reads a version is called from about 1.5–1.6 % (≥ 18–25 reads with a lower bound above 1 %).

The 400× and 1,000× mixtures are subsamples of the same simulated read pools, so they are not independent replicates; the analysis scripts and per-condition results are kept with the validation report.

### Limits

- Haploid model: fractions are shares of genome copies. Polyploid archaea or
  bacteria with several chromosome copies per cell give copy fractions, not
  cell fractions.
- An event inside a repeat longer than the reads (e.g. between two copies of a
  > 20-kb repeat) gets `repeat_limited`.
- Junction estimates weight each molecule type by how often its reads span the
  window. Where versions differ in molecule size (a 40-kb circle versus a
  chromosome), they can differ from the coverage ratio, which weights bases.
- Barcode cross-talk between multiplexed samples places reads of one sample in
  another at roughly 0.1–1 %. Keep `sv_min_fraction` at or above the level
  expected for your kit and run.
- Events below `sv_min_len` (500 bp), and short invertible promoters
  (invertons, typically < 500 bp), are not targeted. Use a dedicated tool
  such as PhaVa for those.
- Discovery depends on Sniffles2 proposing the event. Supplying the genome
  versions removes that dependence for known differences.

## Consensus polishing (s12b_polish)

`s12b_polish` runs `dorado polish --bacteria` on each Autocycler consensus after `s12_collect`,
before any read mapping or annotation uses it. **It is off by default** (`polish_enabled: false`).
In the validation below (rule R1, fixed before any result was seen), polishing this run's drafts at
the depth the rule selects (100×) left the differences to the independent references unchanged or
one higher (131 → 132 summed over 4 genomes), raised the residual errors seen by hold-out reads from
0 to 2 (8 at 50×) and the disagreement between replicate genomes from 7 to 8 sites. The polisher
itself works: it repaired all injected errors (300/300 SNPs and 300/300 homopolymer indels in every
sample at every depth). Autocycler consensus genomes from deep sup reads had nothing left for it to
fix. Turn it on (and re-validate with the scripts in `snv_validation/`) for drafts from shallower
data or a single assembler, or run it in review mode to see where dorado disagrees with a consensus.

### Method

1. **Reads.** A deterministic subsample of `fastq_filtered/<sample>.fastq` at `polish_depth`-fold
   coverage of the draft: a read is kept when a keyed hash of its name falls below the target
   fraction, so the choice does not depend on file order.
2. **Mapping.** `minimap2 -ax lr:hq -y` keeps each read's `RG` tag, and the BAM header gets one `@RG`
   line per read group with `DS:basecall_model=...`; dorado picks its polishing model from it (here
   the bacterial methylation-aware model for sup v5 reads). Reads without `RG` tags are given one
   read group built from `basecall_model`.
3. **Polishing.** `dorado polish --bacteria --any-bam --vcf` on `polish_device` (CPU by default:
   2–5 minutes per 3-Mb genome at 100× with 16 threads on this cluster, most of it one pass over
   the sample's reads).
4. **Change table.** Every record of dorado's VCF (PASS and LowQual: dorado applies both) is listed
   in `changes.tsv`, and applying the table to the draft must reproduce dorado's consensus exactly;
   otherwise nothing is installed. For each change, the table gives the reads whose sequence over the
   change ±5 bp (widened to whole homopolymer runs) matches the draft (`ref_reads`) or the polished
   version (`alt_reads`).
5. **Install.** The polished genome replaces `consensus/<sample>.fasta` with the same contig names and
   FASTA headers (`length=` updated), unless its length differs from the draft by more than
   `polish_max_length_change`. A consensus that s13 mapping or annotation has already used is never
   replaced unless `polish_replace_existing: true`; move those outputs aside and re-run from s12b
   instead. With `polish_install: false` (review mode) everything is written but the consensus is
   kept. A re-run skips samples whose consensus is already the installed polished genome.

### Outputs

`polish/<sample>/`: `draft.fasta` (the Autocycler consensus), `polished.fasta`, `changes.tsv`,
`variants.vcf` (dorado), `info.json` (depth, reads used, change counts, checksums), `stats.json`,
`sub.subsample.json`, `polish.log`, `dorado_polish.log`, `minimap2.log`.

### Validation

All four samples of run `20261012_cm_conjugation` (C. minuta DSM 22607 replicates Cminuta01/02 and
M. smithii T002 conjugation tubes Tube14/15; 2,491–5,074× filtered sup reads) were polished at 50,
100, 200 and 400×. For the drafts and each depth: **Dref**, SNPs + indels ≤ 50 bp to an independent
reference (NCBI NZ_CP029256; the user's T002_i2); **Drep**, the same between the two replicates;
**H**, Clair3 PASS calls from 100× hold-out reads (disjoint from every polishing subsample) mapped to
the genome; **E**, repair of 300 injected SNPs and 300 homopolymer ±1 indels. Rule R1 picks the
smallest depth within 1 of the minimum ΣDref and within 2 of the minimum ΣH, and enables polishing
only if, at that depth, ΣDref, ΣDrep and ΣH are all no worse than for the drafts.

| consensus | ΣDref | ΣDrep | ΣH | injected errors repaired | new differences after repair |
|---|---|---|---|---|---|
| draft | 131 | 7 | 0 | — | — |
| polished 50× | 131 | 7 | 8 | 2,400/2,400 | 0 |
| polished 100× | 132 | 8 | 2 | 2,400/2,400 | 0 |
| polished 200× | 139 | 11 | 1 | 2,400/2,400 | 8 |
| polished 400× | 139 | 13 | 1 | 2,400/2,400 | 10 |

Rule R1 gives `polish_depth: 100` and `polish_enabled: false`; the M. smithii pair alone fails as
well (Dref 14 vs 13, Drep 8 vs 7). Polishing changes fell into two groups: clusters of low-quality
changes that no read in the subsample supports (Cminuta02 at 50–100×: 22 and 16 changes within
1,992,128–1,992,236, QUAL 1–8; Cminuta01 at 200–400×: 4–6 changes near 1,164,700; 117–280 reads
match the draft over each window, 0–2 the polished sequence), and flips at sites where reads are
split roughly evenly between two versions (e.g. Tube15 566,069 G>A, 55 vs 45 reads). The scripts,
the pre-specified rules (`snv_validation/RULES.md`) and the result tables are in `snv_validation/`.

### Settings

| key | default | meaning |
|---|---|---|
| `polish_enabled` | `false` | run s12b_polish (default false: see Validation) |
| `polish_depth` | `100` | fold coverage of the read subsample used for polishing (rule R1: 100) |
| `polish_device` | `cpu` | `cpu`, or `cuda:all` on a GPU node (the job then requests a GPU) |
| `polish_threads` | `16` | threads for minimap2 and dorado polish (SGE slots) |
| `polish_models_dir` | `""` | folder with the dorado polishing models; empty = the folder holding `basecall_model` |
| `polish_max_length_change` | `0.001` | refuse a polished genome whose length differs from the draft by more than this fraction |
| `polish_install` | `true` | `false` = review mode: write polished genome and change table, keep the consensus |
| `polish_replace_existing` | `false` | allow replacing a consensus that s13 mapping or annotation already used |
| `sge_h_vmem_polish` | `6G` | memory per slot |
| `conda_env_snv` | `snv` | env for s12b, s18c and s18d (python ≥ 3.9, mappy, pysam, numpy, scipy, minimap2, samtools, bcftools) |

### Limits

- One run of two species with very deep reads and Autocycler consensus genomes; drafts from
  shallower data or one assembler may carry errors that polishing fixes. Re-validate before enabling.
- NZ_CP029256 differs from both C. minuta consensus genomes at 59 sites whatever the polishing, and
  whether T002_i2 was assembled independently of this run is not known; Dref measures change, not
  absolute accuracy. H and Drep alone already rule out polishing at the selected depth.
- Polishing ignores circularity: reads crossing the start of a circular contig are split, so a few
  hundred bases at contig ends get slightly lower depth.

## Small variants between samples and references (s18c_snv)

`s18c_snv` finds SNPs and short indels (i) between each sample's reads and its own consensus (`self`:
residual consensus errors and majority/minority swaps), (ii) between samples of the run whose
consensus genomes are near-identical (`cross`: replicates, a parent and its mutants) and (iii)
against external genomes listed in `snv_references` (`reference`). The caller configuration and read depth were chosen by validation (rule R2, fixed before any result): **INT** (Clair3 ∩ dorado polish) at **25×** found 3,599/3,600 variants introduced into the tuning genomes and 3,600/3,600 in the held-out genomes, with 0 false positives.

### Method

1. **Plan.** Every pair of consensus genomes (and every sample against every reference) is first
   compared by k-mer containment; plausible pairs are aligned (minimap2 asm5) to get the identity
   over small differences (mismatches and indels ≤ 50 bp; larger indels count as unaligned) and the
   aligned fraction of both genomes. Samples with identity ≥ `snv_cross_min_identity` (0.999) over
   ≥ `snv_cross_min_aligned` (0.9) of both genomes form groups; within a group, the first sample in
   samplesheet order is the reference for the others. A reference genome is used for a sample whose
   consensus matches it at ≥ `snv_ref_min_identity` (0.99) over ≥ `snv_ref_min_aligned` (0.8) of the
   sample. All measured values, used or not, are in `snv/similarity.tsv`.
2. **Reads.** For each comparison, a deterministic `snv_depth`-fold subsample (relative to the
   comparison genome) of the sample's filtered reads, mapped with `minimap2 -ax lr:hq` keeping the
   read groups (as in s12b).
3. **Calling.** The configuration in `snv_caller`: `C3P` = Clair3 `--haploid_precise`; `C3AF` =
   Clair3 diploid model, a heterozygous call becomes the ALT allele if its allele fraction is
   ≥ `snv_clair3_min_af` (0.65) and is dropped otherwise; `DP` = `dorado polish --bacteria --vcf`;
   `INT` / `UNI` = intersection / union of C3AF and DP. PASS calls only, normalised with
   `bcftools norm -m -any`.
4. **Annotation.** For every call: reads supporting the reference and the alternative allele on each
   strand (exact match of the read over the call ±5 bp, widened to whole homopolymer runs), whether it
   lies in a repeat (self-alignment ≥ 300 bp at ≥ 85 % identity), the homopolymer length, the gene and
   its effect (translation table 11; `self`/`cross`: the reference sample's Bakta/Prokka annotation;
   `reference`: the reference's own GFF if given, otherwise the sample's annotation lifted over by
   whole-genome alignment), and `assembly_diff`: whether the two assemblies themselves differ at that
   position. A call where the assemblies agree, or an assembly difference without a call, marks a
   site to inspect.

### Outputs

`snv/comparisons.tsv` (one row per comparison: kind, depth, calls by type, protein-changing calls,
agreement with assembly differences, identity), `snv/calls.tsv` (all calls, columns as above),
`snv/plan.tsv`, `snv/similarity.tsv`, and per comparison `snv/<comparison>/` with `calls.vcf`,
`calls.tsv`, `assembly_diffs.tsv`, `aln.bam` (the subsampled reads; open in IGV with `ref.fasta`),
`info.json` and logs.

Every call carries the gene it falls in, or for an intergenic call the nearest gene and its distance
(`locus_tag`, `gene`, `product`, `effect`, `aa_change`). The report lists every call and a
**Genes with calls** table: one row per comparison and gene with the number of calls, their effects
and the protein changes.

### Validation

Truth came from mutated references: each sample's genome received 1,500 SNPs and 300 indels
(1–20 bp; half change a homopolymer by ±1–2 bp; 10 % in repeat-masked regions), and the sample's own
reads were mapped to the mutated genome at 10, 25, 50 and 100-fold, so every introduced change is a
variant carried by all reads. Five configurations were scored with vcfdist (PASS calls; SNPs and
indels together). Tuning samples Cminuta01/02 chose the configuration; held-out samples Tube14/15
were never used for a choice. Rule R2: highest pooled tuning F1 at 50×, ties (|ΔF1| < 0.0005) to
higher precision, then to a single caller; `snv_depth` = smallest depth within 0.001 of the 100× F1.

| configuration | tuning F1 10× | 25× | 50× | 100× | held-out F1 10× | 25× | 50× | 100× | tuning FP at 50× |
|---|---|---|---|---|---|---|---|---|---|
| C3P (Clair3 haploid_precise) | 0.9978 | 1.0000 | 0.9999 | 0.9997 | 0.9944 | 1.0000 | 1.0000 | 1.0000 | 1 |
| C3AF (Clair3, AF ≥ 0.65) | 0.9985 | 1.0000 | 0.9999 | 0.9997 | 0.9943 | 1.0000 | 1.0000 | 1.0000 | 1 |
| DP (dorado polish) | 0.9602 | 0.9947 | 0.9935 | 0.9976 | 0.9727 | 0.9997 | 0.9997 | 0.9997 | 46 |
| **INT (C3AF ∩ DP), chosen** | 0.9979 | 0.9999 | 0.9999 | 1.0000 | 0.9946 | 1.0000 | 1.0000 | 1.0000 | 0 |
| UNI (C3AF ∪ DP) | 0.9610 | 0.9949 | 0.9935 | 0.9974 | 0.9724 | 0.9997 | 0.9997 | 0.9997 | 47 |

C3P, C3AF and INT tie at 50× (F1 0.99986); INT has the highest precision, so R2 selects
`snv_caller: INT` and `snv_depth: 25`. dorado polish alone makes false positives at every depth
(37 in the tuning samples at 25×); they cluster in the two C. minuta regions where
polishing also changed the drafts without read support (1,150–1,200 kb of Cminuta01, which also
carries many minority variants, and 1,950–2,000 kb of Cminuta02). Requiring both callers removes them.
At 25× and above, INT found every introduced SNP and indel inside repeat-masked regions. Below 25×,
every configuration misses variants. Scripts, rules and tables: `snv_validation/`.

### Settings

| key | default | meaning |
|---|---|---|
| `snv_enabled` | `true` | run s18c_snv |
| `snv_caller` | `INT` | `C3P`, `C3AF`, `DP`, `INT` or `UNI` (rule R2: `INT`) |
| `snv_depth` | `25` | fold coverage of the read subsample per comparison (rule R2: 25) |
| `snv_self` | `true` | compare each sample with its own consensus |
| `snv_cross` | `true` | compare samples whose consensus genomes are near-identical |
| `snv_cross_min_identity` | `0.999` | identity over small differences needed to group samples |
| `snv_cross_min_aligned` | `0.9` | aligned fraction of both genomes needed to group samples |
| `snv_references` | `[]` | external genomes: FASTA paths or `{fasta, gff, name}` entries |
| `snv_ref_min_identity` | `0.99` | identity needed to use a reference for a sample |
| `snv_ref_min_aligned` | `0.8` | aligned fraction of the sample needed to use a reference |
| `snv_clair3_min_af` | `0.65` | C3AF/INT/UNI: allele fraction at which a Clair3 heterozygous call becomes the ALT allele |
| `conda_env_clair3` | `clair3` | env name or absolute prefix holding Clair3 ≥ 2.0 |
| `snv_clair3_model` | `/tmp/global2/jmarsh/code/conda_envs/clair3_models/r1041_e82_400bps_sup_v520` | Clair3 model folder (r1041_e82_400bps_sup_v520) |
| `snv_threads` | `16` | threads per comparison (SGE slots) |
| `sge_h_vmem_snv` | `6G` | memory per slot |

### Limits

- Haploid calling: a site carried by part of the population is reported only once it is the
  majority in the subsample (C3P) or above `snv_clair3_min_af` (C3AF). Minority variants are the job
  of s18d_minor.
- Repeats longer than the reads (rRNA operons, IS copies) are flagged (`in_repeat`), not resolved.
- Liftover covers genes in alignable regions only; genes absent from the sample are not annotated on
  the reference.

## Minority variants within samples (s18d_minor)

`s18d_minor` looks for single-nucleotide variants carried by part of a sample's population, using
the reads of the s13 mapping. Each sample is tested against its own error profile, learned from the
same reads, so systematic basecalling errors of a sequence context are absorbed into the error model
instead of being reported as variants.
Defaults were fixed before validation and accepted by rule R3 in all 48 spike-in conditions (two samples, 100–2,000×, 0.25–40 %): at most 6 false calls per condition, the 95 % interval contained the true fraction for 6,941 of 6,942 called spiked sites (true fraction ≥ 1 %), and no site inside a masked repeat was called. At ~2,000× most minority variants are found from 2 %, at ~500× from 5 %.

### Method

1. **Counts.** Bases per position and strand from reads with MAPQ ≥ `minor_min_mapq` (20) and base
   quality ≥ `minor_min_baseq` (10). Positions in repeats (self-alignment ≥ 300 bp at ≥ 85 %
   identity), within `minor_end_mask` of a linear contig end, or with fewer than
   `minor_min_strand_depth` reads on either strand are not assessed.
2. **Error model.** For each strand, alternative base and sequence context (7-mer, 5-mer, 3-mer and
   1-mer, each combined with the homopolymer length), a beta-binomial error rate is fitted to all
   assessed positions of the sample after trimming outlying positions; the most specific context
   with enough positions is used.
3. **Test.** A position-base pair is a candidate when both strands are significant (p = the larger of
   the two one-sided beta-binomial p-values) and each strand carries ≥ `minor_min_alt_per_strand`
   ALT reads. Benjamini-Hochberg over all positions × 3 bases, q ≤ `minor_fdr`. Candidates with strong
   strand imbalance, or in a context enriched in candidates genome-wide (a systematic error motif),
   are removed.
4. **Fraction.** Error-corrected fraction with a Wilson 95 % interval; the detection floor is the
   smallest fraction that would reach significance at a typical assessed position.
5. **Linkage.** Sites carried by the same reads (odds ratio ≥ 20, Fisher p ≤ 1e-8, ≤ 50 kb apart)
   form groups; weaker sites on the reads of a group are reported as `linked`. A group's fraction and
   interval come from the reads covering its sites.
6. **Known versions.** If `sv_variant_genomes` lists two or more genome versions, their SNP
   differences are genotyped in every sample (`known_sites.tsv`, `known_versions.tsv`).

### Outputs

`minor/<sample>/sites.tsv` (one row per site: fraction and 95 % interval, raw fraction, background
error rate, ALT reads and depth per strand, q, context, homopolymer length, group, gene and effect),
`groups.tsv` (one row per group of linked sites, last column `genes`: span, sites, reads, fraction
and interval, positions, and the genes of the group's sites as `locus_tag(gene):product [class xN, …]`
joined by ` | `, where class is the site effect or `intergenic`, an intergenic site counting for its
nearest gene), `summary.tsv` (positions assessed and masked, median depth, detection floor, counts),
`null_model.tsv` (error rates per context level), and `minor/summary.tsv` for all samples.

The report lists every site with its gene and effect, every group with its genes, and a
**Genes with minority variants** table: one row per sample and gene with the number of sites (called
and linked), their effects, the highest fraction and each protein change with its fraction
(intergenic sites count for the nearest gene).

### Validation

Real reads of Cminuta01 and Tube14 (mapped as in s13) were subsampled to 100, 500 and 2,000×. A
minor haplotype (a random set of reads, 0.25–40 % of all reads) received a non-reference base at 300
sites (200 isolated, 20 clusters of 5) plus 30 sites in masked repeats; each condition and the
unmodified reads were analysed with the production defaults. Rule R3 (RULES.md): false calls
(non-spiked positions not called in the unmodified reads) ≤ 5 % of calls (≤ 1 when fewer than 20
calls), ≥ 90 % interval coverage for spiked sites ≥ 1 %, no masked site called. Spiked sites
reported (called or linked):

| sample | depth | 0.25 % | 0.5 % | 1 % | 2 % | 5 % | 10 % | 20 % | 40 % | false calls (max per condition) | background calls in control |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C. minuta 01 | 100× (median 96×) | 0/300 | 0/300 | 0/300 | 0/300 | 3/300 | 120/300 | 287/300 | 300/300 | 1 | 101 |
| C. minuta 01 | 500× (median 462×) | 0/300 | 2/300 | 3/300 | 54/300 | 268/300 | 296/300 | 300/300 | 300/300 | 0 | 277 |
| C. minuta 01 | 2,000× (median 1,848×) | 0/300 | 3/300 | 103/300 | 255/300 | 296/300 | 298/300 | 300/300 | 300/300 | 6 | 374 |
| M. smithii Tube14 | 100× (median 94×) | 0/300 | 0/300 | 0/300 | 0/300 | 6/300 | 89/300 | 283/300 | 299/300 | 0 | 1 |
| M. smithii Tube14 | 500× (median 485×) | 0/300 | 0/300 | 1/300 | 48/300 | 277/300 | 299/300 | 299/300 | 299/300 | 0 | 3 |
| M. smithii Tube14 | 2,000× (median 1,945×) | 0/300 | 21/300 | 183/300 | 280/300 | 297/300 | 299/300 | 299/300 | 299/300 | 1 | 10 |

R3 holds in all 48 conditions: 42 false calls in total (≤ 6 per condition,
false-discovery proportion ≤ 3.1 %), interval coverage 6,941/6,942, masked spiked sites
called 0/1,440, fraction error ≤ 0.22 percentage points. The unmodified reads already
carry minority variants (Cminuta01: 101–374 sites in seven linked clusters; Tube14: 1–10); these are
reported, not counted as false. Production runs subsample deeper samples to `minor_max_depth`
(2,000×), the top of the validated range. Scripts, rules and tables: `snv_validation/`.

### Settings

| key | default | meaning |
|---|---|---|
| `minor_enabled` | `true` | run s18d_minor (rule R3 met: true) |
| `minor_fdr` | `0.01` | Benjamini-Hochberg q threshold |
| `minor_min_alt_per_strand` | `3` | ALT reads needed on each strand |
| `minor_min_strand_depth` | `20` | reads needed on each strand for a position to be assessed |
| `minor_min_baseq` | `10` | minimum base quality counted |
| `minor_min_mapq` | `20` | minimum mapping quality counted |
| `minor_end_mask` | `500` | bases masked at linear contig ends |
| `minor_max_depth` | `2000` | deeper samples are read-subsampled to this mean depth (top of the validated range; 0 = all reads) |
| `minor_parallel` | `4` | samples analysed at the same time (one SGE slot each) |
| `sge_h_vmem_minor` | `8G` | memory per slot |

### Limits

- SNVs only; minority indels are dominated by homopolymer errors in nanopore reads and are not
  called.
- Barcode cross-talk places 0.1–1 % of another sample's reads in each sample: a minority variant
  near that level that matches another sample of the run should be treated as possible cross-talk.
- Repeats and contig ends are masked, so variants there are not assessed.

## Phase variation (s18e_phava)

`s18e_phava` finds invertible DNA segments (invertons): stretches between inverted repeats that part
of a population carries in the opposite orientation, such as invertible promoters that switch genes
on and off (phase variation). It measures the inverted fraction in every sample with PhaVa 0.2.3
(Chanin, West et al. 2024, *Nature* 634: 234–242, doi:10.1038/s41586-024-07970-4) and gives each inverton the genes it lies in, cuts or
could drive. The stage is optional (`phava_enabled`; off when the key is absent, on in the shipped
config).
Defaults were fixed before validation and accepted by rule R4 in both samples: at 5, 20 and 50 % inverted reads all 20 spiked invertons of each sample were called, every fraction was within 0.05 of the truth (median error at most 0.008), and the unspiked invertons called were exactly those of the control. PhaVa's own joint counts overestimated 6 of the 40 spiked invertons by 0.06–0.43, which is why every candidate with an inverted read is re-counted (step 5).

### Method

1. **Annotation for PhaVa.** `phava_tools.py gff` turns the sample's Bakta or Prokka GFF into the
   GFF that PhaVa 0.2.3 reads (`phava_genes.gff`): a `region` line per consensus contig, every
   non-pseudo CDS on its consensus contig name (contigs matched by sequence, so Bakta's renamed
   contigs are found), `locus_tag` and `transl_table` (11 when absent) on every CDS, and no FASTA
   section. PhaVa stops with an error on the unconverted Bakta and Prokka files.
2. **Reads.** The filtered reads (`fastq_filtered/<sample>.fastq`) are subsampled by a keyed hash
   of the read name to `phava_max_depth` when the sample is deeper (`subsample.json`; the method of
   s12b, s18c and s18d). PhaVa's alignment file and run time grow with depth.
3. **Candidates.** PhaVa runs EMBOSS einverted twice on the consensus (thresholds 51 and 75, repeat
   arms up to 750 bp, gaps up to 100 bp) and keeps the inverted repeats with at least 30 bp between
   them. For every candidate it builds the reference and the inverted orientation (the segment
   between the repeats reverse-complemented) with 1 kb of flank on each side.
4. **Counting.** minimap2 maps the reads to both orientations of all candidates. A read counts for
   an orientation when it aligns there with MAPQ above 2 and at most 15 % mismatches in the segment
   and in the flanks, does not also align to the other orientation, and spans the whole inverton
   with more than 20 bp on each side. `variation_wf --reportAll` reports every candidate.
5. **Re-count** (`phava_tools.py genotype`). PhaVa maps the reads to the references of all
   candidates at once. Where candidates are close, their references share flank sequence: a read
   in the reference orientation then fits two references equally well (MAPQ 0) and is not
   counted, while a read with the inverted segment fits only one. The inverted fraction is then
   inflated. In Cminuta01, PhaVa counted 23 of 89 reads inverted at 1:246,611 (25.8 %), against
   24 of 808 (3.0 %) for the same reads counted at that inverton alone; the s18d minority group
   there gives 3.0 %. Inverted reads are not lost this way, so every candidate with at least one
   PhaVa inverted read is re-counted against its own two references only. The re-count uses
   PhaVa's references, minimap2 command and counting rules, and the same reads, which are found
   through the s13 mapping (`genotype.tsv`). Calls and fractions use these counts
   (`counts_from` = `regenotyped`); PhaVa's joint counts are kept (`phava_forward_reads`,
   `phava_reverse_reads`). Without an s13 mapping the stage keeps PhaVa's counts and logs a
   warning.
6. **Call.** An inverton is **called** when at least `phava_min_reverse_reads` reads carry the
   inverted orientation and they make up at least `phava_min_reverse_fraction` of the counted reads
   (PhaVa's recommended filter). The inverted fraction is reverse / (forward + reverse) with a
   Wilson 95 % interval. A candidate that no read spans gets `NA` (PhaVa's own table writes 1 there).
7. **Gene context** (`phava_tools.py annotate`; features are all GFF lines except gene, region and
   source, labelled `locus_tag(gene):product` as in the other stages):
   - `class`: the segment against the CDS, as PhaVa does: `intragenic` (inside one CDS; the
     inversion recodes or truncates it), `gene_boundary` (overlaps a CDS in part), `intergenic`.
   - `genes_in_segment` (every feature overlapping the segment) and `breakpoint_genes` (the
     features at both ends of the segment, in the notation of s18b:
     `start C:S in X [cut] | end C:E intergenic (left A, n bp; right B, m bp)`; `[cut]` = the
     feature extends beyond the segment).
   - `left_gene` and `right_gene`: the nearest gene (CDS, tRNA, rRNA, ncRNA, tmRNA) entirely
     outside the inverton on each side, with strand and the bases between it and the inverton
     (wrapping round circular contigs).
   - `promoter_targets`: for an intergenic inverton, a flanking gene whose 5′ end faces the inverton
     within 300 bp (right gene on +, left gene on −), i.e. a gene an invertible promoter in the
     segment could switch. This is position only, not evidence of a promoter.
8. **Across samples.** `sequence_key` identifies the segment sequence in either orientation, so
   the same inverton is linked between samples even when one assembly is rotated or reverse
   complemented (`in_other_samples` of `invertons_called.tsv`).

### Outputs and statuses

Per sample, `phava/<sample>/`:

- `invertons.tsv`: one row per PhaVa candidate, in contig order.
- `genotype.tsv`: the re-count of the candidates with PhaVa inverted reads (`reads_tested`, re-counted
  and PhaVa forward and reverse reads).
- `summary.tsv`: counts for the sample.
- `phava_ratio.tsv`: PhaVa's own table.
- `phava_IRs.tsv`: the candidates with their repeat and segment sequences.
- `phava_genes.gff`: the GFF PhaVa used.
- `phava_genes.json`: the conversion counts.
- `subsample.json`: the reads used.
- `phava.log`: PhaVa's log.

Run level: `phava/summary.tsv` (all samples), `phava/invertons_called.tsv` (the called invertons of
all samples, plus `in_other_samples`) and `phava/invertons_all.tsv` (every candidate of every sample).

| column | meaning |
|---|---|
| `sample` | sample |
| `inverton` | PhaVa identifier `contig:a-b-c-d` (0-based starts, 1-based ends: left repeat a+1..b, segment b+1..c, right repeat c+1..d) |
| `contig`, `start`, `end` | the inverton, repeats included (1-based, inclusive) |
| `left_ir`, `right_ir`, `segment` | the two inverted repeats and the invertible segment between them |
| `segment_bp`, `ir_bp` | lengths of the segment and of the left repeat |
| `forward_reads`, `reverse_reads`, `counted_reads` | reads PhaVa counted in the reference and the inverted orientation, and both |
| `reverse_fraction`, `ci_low`, `ci_high` | inverted fraction and Wilson 95 % interval (`NA` without counted reads) |
| `status` | `called`; `below_filter` (inverted reads, but fewer than the thresholds); `no_reverse_reads`; `not_spanned` (no read counted) |
| `counts_from` | `regenotyped` (re-counted alone, step 5) or `phava` (PhaVa's joint count; candidates without a PhaVa inverted read) |
| `phava_forward_reads`, `phava_reverse_reads` | PhaVa's own joint counts |
| `class` | `intragenic`, `gene_boundary` or `intergenic` (segment against the CDS) |
| `genes_in_segment`, `breakpoint_genes` | features overlapping the segment; features at its two ends |
| `left_gene`, `left_strand`, `left_distance` | nearest gene left of the inverton, its strand and the bases in between |
| `right_gene`, `right_strand`, `right_distance` | the same on the right |
| `promoter_targets` | flanking genes whose 5′ end faces an intergenic inverton within 300 bp, with strand and distance |
| `phava_gene_overlaps`, `phava_upstream`, `phava_downstream` | PhaVa's own gene columns (effect class of intragenic invertons; nearest genes by coordinate) |
| `sequence_key`, `key_copies` | segment sequence identifier (either orientation) and the number of candidates of the sample with it |
| `overlapping_candidates` | other candidates of the sample whose span overlaps this one (one inversion can be reported at more than one candidate) |
| `in_other_samples` | (`invertons_called.tsv`) the same segment in the other samples: status, fraction and reads |

`summary.tsv` columns: `sample`, `candidates`, `spanned` (at least one counted read),
`with_reverse_reads`, `regenotyped` (candidates re-counted), `phava_called` (calls PhaVa's joint
counts would give), `called`, `called_intergenic`, `called_gene_boundary`, `called_intragenic`,
`called_with_promoter_target`, `median_counted_reads`, `reads_used`, `depth_used`, `subsampled`,
`min_reverse_reads`, `min_reverse_fraction`, `annotation` (the GFF used or `none`),
`phava_version` and `tools_version`.

The report's **Phase variation** section has these tables:

- the per-sample summary;
- **Genes at called invertons**: one row per called inverton and gene (genes in the segment,
  promoter targets, and for other intergenic invertons the nearest gene on each side);
- **Called invertons**, with PhaVa's joint counts beside the re-count;
- **Inverted reads below the call threshold**;
- with s18d enabled, **Minority-variant groups and candidate invertons**: every s18d group with the
  candidate holding most of its sites.

Reads carrying an inverted segment differ from the consensus at many linked positions inside it, so
an s18d group inside a candidate inverton may be that inversion rather than point mutations.

### Validation

Real reads of Cminuta01 and Tube14 (the production subsample at 1,000×) were analysed unmodified
(control) and with 20 invertons per sample spiked. In 5, 20 or 50 % of the reads whose s13 alignment spans an
inverton with 50 bp to spare, the bases aligned to the segment were reverse-complemented and their qualities
reversed, which is the inverted orientation PhaVa tests. The spiked invertons were drawn (seed 1) from the
control candidates with no inverted read, at least 100 counted reads and a 30–2,000 bp segment, overlapping
no other candidate, and at least 10 kb from each other and from any candidate with inverted reads. Each read
set was analysed like production, with PhaVa in its own directory. Rule R4 (`phava_validation/RULES.md`):

- at least 19 of 20 spiked invertons are called in every condition;
- the inverted fraction is within 0.05 of the truth for at least 90 % of them, with a median error of at
  most 0.02;
- the called unspiked invertons equal the control's;
- no spiked inverton has inverted reads in the control.

| sample | inverted reads spiked | called | within 0.05 | median error | largest error | PhaVa joint count within 0.05 | other invertons called (control) |
|---|---|---|---|---|---|---|---|
| C. minuta 01 | 5% | 20/20 | 20/20 | 0.0004 | 0.0020 | 19/20 | 3 (3) |
| C. minuta 01 | 20% | 20/20 | 20/20 | 0.0023 | 0.0300 | 19/20 | 3 (3) |
| C. minuta 01 | 50% | 20/20 | 20/20 | 0.0054 | 0.0301 | 19/20 | 3 (3) |
| M. smithii Tube14 | 5% | 20/20 | 20/20 | 0.0007 | 0.0035 | 15/20 | 0 (0) |
| M. smithii Tube14 | 20% | 20/20 | 20/20 | 0.0035 | 0.0064 | 15/20 | 0 (0) |
| M. smithii Tube14 | 50% | 20/20 | 20/20 | 0.0077 | 0.0190 | 15/20 | 0 (0) |

R4 holds in all six conditions. PhaVa's joint counts alone would have failed the accuracy criterion in
Tube14 (15 of 20 within 0.05). At the 6 invertons it
overestimated (1 in C. minuta 01, 5 in M. smithii Tube14), PhaVa counted only 133–662 reads, against a
median of 908 in the re-count: reads in the reference orientation were lost in its joint mapping. At every
re-counted inverton of the run, the re-count also agrees to within 0.0014 with an independent count of the
same reads (mappy, one mismatch rule over the alignment). The stage takes about 1 h per sample at 1,000×
with 4 threads (2.2 GB memory, 9–11 GB temporary alignment file per sample). Scripts, rules and tables are
in `phava_validation/`.

### Settings

| key | default | meaning |
|---|---|---|
| `phava_enabled` | `true` | run s18e_phava (off when the key is absent) |
| `conda_env_phava` | `phava` | conda env with phava 0.2.3 (bioconda; brings EMBOSS einverted and minimap2) |
| `phava_max_depth` | `1000` | deeper samples are read-subsampled to this mean depth (0 = all reads) |
| `phava_min_reverse_reads` | `3` | inverted reads needed for a call |
| `phava_min_reverse_fraction` | `0.03` | inverted / counted reads needed for a call |
| `phava_threads` | `4` | threads per sample |
| `phava_parallel` | `4` | samples analysed at the same time (the stage requests threads × parallel slots) |
| `sge_h_vmem_phava` | `4G` | memory per slot |

### Limits

- einverted searches the linear sequence, so an inverton across the origin of a circular contig is
  not found, and only repeats einverted finds are candidates (arms up to 750 bp, segments of at
  least 30 bp).
- A read must span the whole inverton plus 20 bp on each side, so segments approaching the read
  length get few or no counted reads (`not_spanned`).
- An inverted fraction below `phava_min_reverse_fraction` (3 %) is not called, whatever the depth;
  barcode cross-talk (0.1–1 % of another sample's reads) stays below it.
- PhaVa keeps overlapping candidates (in *M. smithii* about one in eight), so one inversion can be
  called at two candidates; `overlapping_candidates` names them.

## Assumptions

- The conda envs named in `pipeline_config.yaml` exist and contain the
  expected tools: `nanopore` → dorado + chopper + flye, `autocycler` →
  autocycler plus every assembler in `assemblers:`, `prokka` → Bakta,
  `gtdbtk3` → GTDB-Tk (with `GTDBTK_DATA_PATH` set), `mmseqs` → minimap2 +
  samtools + metabat2 + pandas/matplotlib, `checkm2` → CheckM2, and (if
  using MIJAMP) `mijamp` → pandas, meme, ont-modkit, minimap2, samtools,
  biopython, xmltodict, and `sv` → python ≥ 3.9, sniffles ≥ 2.2, mappy,
  pysam, samtools; `snv` → python ≥ 3.9, mappy, pysam, numpy, scipy,
  minimap2, samtools, bcftools (vcfdist only for re-running the validation);
  `clair3` → Clair3 ≥ 2.0 with the model in `snv_clair3_model`.

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

`s21_mijamp` runs MIJAMP's **preprocessing** step only (`scripts/preprocess`:
minimap2 alignment of each sample's demultiplexed BAM to its genome, then
modkit), leaving `mijamp/<sample>/genome.fasta`, `mapped.bam`,
`mapped.bam.bai` and `data.bed`.

Two settings in `pipeline_config.yaml` prepare the data for MIJAMP's `motif`
step:

- `mijamp_target_depth` (default `100`; `0` keeps all reads): each sample's
  BAM is downsampled to about this read depth before preprocessing.
- `mijamp_filter_threshold` (default `0.9`): passed to modkit as
  `--filter-threshold` (MIJAMP's own default is 0.6).

`motif` sends a site to MEME only if it is ≥ 90 % modified and no read shows
a different base there (`Ndiff` = 0), and it needs more than 100 such sites
per modification type. These filters were set for ~100× data (MIJAMP's test
data are 100×). At several hundred × and above nearly every site has a few
reads with sequencing errors, and with a threshold of 0.6 many fully
methylated 6mA sites stay below 90 %, so `motif` reports no motifs even when
the genome is methylated.

Motif discovery is **not** run by the pipeline, and `motif` is interactive
(it asks you to confirm each motif). Run it per sample when you need it:

```bash
mamba activate mijamp
<mijamp_dir>/scripts/motif -f <output-dir>/mijamp/<sample>/ -t 20
```

MIJAMP writes the motif table to `output.tsv` in the same folder.

- Dorado demux output filenames contain the barcode label. The `s03_rename`
  step matches `*<barcode>*.fastq[.gz]`.
- Autocycler subsample produces 4 subsamples named `sample_01.fastq` …
  `sample_04.fastq`.
- SGE job names need a leading letter — stages are named `s01_*` … `s13_*`
  to satisfy this.
