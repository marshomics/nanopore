# Nanopore assembly + annotation pipeline

End-to-end SGE pipeline that takes raw POD5 files and a sample sheet, runs
basecalling → demux → QC → Autocycler 8-assembler consensus → Prokka, and
ends with one annotated assembly per sample.

## Files

- `generate_pipeline.py` — emits all SGE submit scripts for a given run
- `pipeline_config.yaml` — paths, models, conda envs, SGE resources, filter thresholds
- `samplesheet_template.tsv` — fill in one row per barcode

## Sample sheet

TSV with these columns (header required, order doesn't matter):

| column              | meaning                                              |
|---------------------|------------------------------------------------------|
| `barcode`           | dorado demux barcode tag, e.g. `barcode01`           |
| `sample`            | output prefix; must be unique                        |
| `genome_length`     | expected genome size in bp (integer; commas allowed) |
| `genus`             | Prokka `--genus`                                     |
| `species`           | Prokka `--species`                                   |
| `strain`            | Prokka `--strain` (and `--prefix` / `--locustag`)    |
| `kingdom`           | `Bacteria` or `Archaea`                              |
| `reference_proteins`| path to a `.gbk`; leave blank to skip `--proteins`   |

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
| `s13_annotate`    | `prokka/<sample>/<sample>.gff` exists                                           |

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
├── prokka/<sample>/               # s13_annotate
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
                                              └─ s13_annotate (Prokka)
```

## Assumptions

- The conda envs named in `pipeline_config.yaml` exist and contain the
  expected tools: `nanopore` → dorado + chopper, `autocycler` → autocycler
  plus every assembler in `assemblers:`, `prokka` → prokka.
- Dorado demux output filenames contain the barcode label. The `s03_rename`
  step matches `*<barcode>*.fastq[.gz]`.
- Autocycler subsample produces 4 subsamples named `sample_01.fastq` …
  `sample_04.fastq`.
- SGE job names need a leading letter — stages are named `s01_*` … `s13_*`
  to satisfy this.
