# Validation of s12b_polish, s18c_snv and s18d_minor

The defaults of the three small-variant stages were set by rules fixed in
`RULES.md` before any result was seen. Every rule was applied to the tables in
`results/`, and the gate checks that recomputed those numbers from the raw tool
outputs are kept with the session notes. Data: run `20261012_cm_conjugation`
(Cminuta01/02: *Christensenella minuta* DSM 22607; Tube14/15_T002_conj:
*Methanobrevibacter smithii* T002), sup v5.2.0 reads, dorado 2.1.2, Clair3 2.0.3,
vcfdist 2.6.4, minimap2 2.31, samtools/bcftools 1.24.

| file | what it does |
|---|---|
| `install_envs.sh` | creates the `clair3` and `snv` envs and downloads the Clair3 model; records versions and checksums |
| `prep_refs.py` | writes NZ_CP029256 (from the run's genomic.gbff) and T002_i2 as FASTA references |
| `polval.sge` | rule R1: polishes each sample at 50/100/200/400×, hold-out Clair3 calls, differences to the references, injected-error repair (SGE array, one task per sample) |
| `snvval.sge` | rule R2: mutated reference per sample, real reads at 10/25/50/100×, five caller configurations, vcfdist (`qsub -v TDEPTH=100`) |
| `minprep.sge`, `minval.sge` | rule R3: all reads mapped as in s13, spike-in sites; then the control and 8 spike-in fractions at 100/500/2,000× per sample |
| `minval_one.sge` | one spike-in fraction of one condition, same commands and seeds as `minval.sge` (used to run the Cminuta01 2,000× fractions in parallel) |
| `spike_minor.py` | chooses spike-in sites and edits the minor reads' bases in a BAM |
| `val_aggregate.py`, `agg.sge` | collects each validation into the tables in `results/` |

Paths at the top of each script point to the cluster locations used here; change
them to re-run on other data. The engines themselves (`polish_tools.py`,
`snv_tools.py`, `minor_variants.py`) are the ones next to `generate_pipeline.py`.

Results: `results/polishing/` (summary.tsv: Dref and H per sample and state;
pairs.tsv: Drep; injections.tsv), `results/small_variants/` (summary.tsv: vcfdist
counts and F1 per sample, depth and configuration; strata.tsv: recall by variant
kind and repeat status), `results/minority_variants/` (conditions.tsv, summary.tsv:
calls, false calls, background and interval coverage per condition;
sensitivity.tsv: called and linked spike-in sites per set). Figures: one PDF per
stage in the same folders.
