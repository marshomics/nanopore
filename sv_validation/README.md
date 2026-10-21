# s18b_sv validation harness

Scripts and results behind the "Validation" table in the main README. Paths
inside the SGE scripts point at the scratch area where the validation ran
(`/ebio/abt3_scratch/jmarsh/software/science/sv_build`); edit `S=` / the
`#$ -o` lines to rerun elsewhere. Everything runs in the `sv` conda env
(Badread 0.4.2, minimap2 2.31, samtools 1.24, Sniffles2 2.8.1, mappy/pysam).

| file | role |
|---|---|
| `make_synthetic_versions.py` | builds `T002_i2_inv3k` (inversion of 400,001–403,000) and `T002_i2_dup5k` (tandem duplication of 1,400,001–1,405,000) from `T002_i2.fasta` |
| `truth_events.json` | the three events, in base-genome coordinates |
| `sim_{A,B,INV,DUP}.sh` | Badread read pools: `--quantity 50x` per array task (A, B: 20 tasks = 1,000×; INV, DUP: 8 tasks = 400×), `--length 11000,10000 --identity 98,99.8,2 --error_model nanopore2023 --qscore_model nanopore2023`; read names are prefixed with the source version |
| `map_<pool>_on_<consensus>.sh` | `minimap2 -ax map-ont` of each pool to the consensus it is tested on (A = `T002_i2`, B = `T002_i2_with_deletion`); also records each pool's total bases |
| `conditions.tsv` | the 120 conditions (event, consensus, depth, fraction) |
| `val_array.sh` | per condition: subsample the two pool BAMs to the target depth and fraction (`samtools view -s`, fixed seeds), merge, run `sv_genotype.py` with both genome versions and again with discovery only, score with `sv_eval.py` |
| `sv_eval.py` | scores one condition against the realised truth (share of reads from each version) and checks every assigned read against its name prefix |
| `analyse_validation.py` | tables + figure from the per-condition JSONs |
| `validation_results.tsv`, `validation_summary.tsv`, `validation_accuracy.pdf` | results for engine version 1.2; version 1.3 differs only for loci with allele windows > `sv_max_span_window` (20 kb), and no validation locus has one (largest window 12 kb), so these results hold for 1.3 |

Mixtures at 400× and 1,000× are subsamples of the same read pools (not
independent replicates), and the consensus-A and consensus-B conditions of
the deletion use the same reads mapped to the two consensuses.
