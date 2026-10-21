# Validation of s18e_phava

The settings of the phase-variation stage (PhaVa 0.2.3, read cap, call thresholds) were tested
under rule R4, fixed in `RULES.md` before any spiked result was seen. The amendments are dated
there: A1 before any PhaVa run, A2 after the two control runs. The gate checks that recomputed
every number from the raw PhaVa tables are kept with the session notes. Data: run
`20261012_cm_conjugation` (Cminuta01: *Christensenella minuta* DSM 22607; Tube14_T002_conj:
*Methanobrevibacter smithii* T002), sup v5.2.0 reads, minimap2 2.31, EMBOSS einverted 6.6.0.0,
pysam 0.24.1.

| file | what it does |
|---|---|
| `time_phava.sge` | timing run at the read cap: GFF conversion, read subsample, `phava variation_wf --reportAll`, annotation (time, memory and disk recorded); its PhaVa run is the R4 control |
| `val_prep_from_timing.sge` | rule R4 per sample: takes the control from the timing run, draws the 20 spiked invertons and writes the edited read sets |
| `val_prep.sge` | the same with its own control run (for re-running from scratch) |
| `val_cond.sge` | PhaVa on one spiked read set, in its own directory |
| `geno_val.sge` | re-count (`phava_tools.py genotype`) and annotation of one condition, as in production |
| `spike_inv.py` | the draw (`select`) and the read edits (`edit`): in a seeded fraction of the reads whose s13 alignment spans an inverton with 50 bp to spare, the bases aligned to the segment are reverse-complemented and their qualities reversed |
| `val_metrics.py` | metrics per sample and fraction (`results/metrics.tsv`) and per spiked inverton (`results/val_invertons.tsv`) |
| `plot_accuracy.py` | the figure `results/phava_validation_accuracy.png` |
| `inv_regeno.py`, `regeno.sge` | independent count of an inverton's two orientations with mappy (cross-check of the re-count) |

Paths at the top of each script point to the cluster locations used here; change them to re-run
on other data. The engine (`phava_tools.py`) is the one next to `generate_pipeline.py`.

Results (`results/`):

- `metrics.tsv`: one row per sample and fraction, giving called spiked invertons, how many are within 0.05,
  the median and largest error, the same for PhaVa's joint count, and the unspiked calls against the
  control.
- `val_invertons.tsv`: every spiked inverton and fraction, with truth, re-counted and PhaVa counts.
- `spiked_*.tsv`, `truth_*.tsv`: the draw and the edited reads.
- `design.json`: the settings, written before the first spiked run.
- `timing.json`: run time, memory and disk at the cap.

Rule R4 holds in all six conditions. All 20 spiked invertons are called in each, every re-counted
fraction is within 0.05 of the truth (median error 0.0004–0.0077), and the unspiked calls equal the
control's (Cminuta01: the same 3; Tube14: none). PhaVa's joint counts were off by more than 0.05 at 6 of
the 40 spiked invertons (overestimates of 0.06–0.43). With them, Tube14 would have failed the accuracy
criterion (15 of 20 within 0.05); this is amendment A2's re-count.

