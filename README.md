# Nanopore assembly + annotation pipeline

End-to-end SGE pipeline that takes raw POD5 files and a sample sheet, runs
basecalling → demux → QC → Autocycler 8-assembler consensus → Prokka, and
ends with one annotated assembly per sample.

## Files

- `generate_pipeline.py` — emits all SGE submit scripts for a given run
- `pipeline_config.yaml` — paths, models, conda envs, SGE resource defaults
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
    --pod5-dir /ebio/abt3_scratch/jmarsh/nanopore/no_sample_id/<run_id>/pod5 \
    --samplesheet /path/to/samplesheet.tsv \
    --output-dir /ebio/abt3_scratch/jmarsh/nanopore/<run_name> \
    --config pipeline_config.yaml

# 4. Submit. Every stage holds on the previous via -hold_jid.
bash <output-dir>/submit_all.sh
```

## Output layout

```
<output-dir>/
├── basecalled/basecalls.bam       # 01_basecall
├── fastq_raw/                     # 02_demux (dorado output)
├── fastq/<sample>.fastq           # 03_rename
├── fastq_filtered/<sample>.fastq  # 04_filter (chopper)
├── autocycler/<sample>/           # 05_subsample, 08_compress, 09_cluster
├── autocycler/assemblies/         # 06_assemble (flat)
├── autocycler/assemblies_fasta/   # 07_organize (per-sample)
├── autocycler/cluster/<sample>/   # 09_cluster (symlink)
├── consensus/<sample>.fasta       # 12_collect
├── prokka/<sample>/               # 13_annotate
├── submit_scripts/                # generated SGE scripts
├── stdout/                        # SGE log files
└── submit_all.sh
```

## Stage order and hold dependencies

```
01_basecall (GPU)
  └─ 02_demux (GPU)
      └─ 03_rename
          └─ 04_filter (chopper)
              └─ 05_subsample
                  └─ 06_assemble_batch_001..NNN  (parallel, jobs_per_assembly_batch)
                      └─ 07_organize
                          └─ 08_compress
                              └─ 09_cluster
                                  └─ 10_trim_resolve
                                      └─ 11_combine
                                          └─ 12_collect
                                              └─ 13_annotate (Prokka)
```

## Editing

If you need to change a single stage, edit the generated `<output-dir>/submit_scripts/NN_*.sh`
and resubmit just that script. Hold-job names follow the script filenames
(without `.sh`), so you can re-wire dependencies manually with
`qsub -hold_jid <name> ...`.

## Assumptions

- The conda envs named in `pipeline_config.yaml` exist and contain the
  expected tools: `nanopore` → dorado + chopper, `autocycler` → autocycler
  + every assembler in `assemblers:`, `prokka` → prokka.
- Dorado demux output filenames contain the barcode label (default since
  dorado 0.5+). The 03_rename step matches `*<barcode>*.fastq[.gz]`.
- Autocycler subsample produces 4 subsamples named `sample_01.fastq` …
  `sample_04.fastq` — the assembly stage assumes this.
