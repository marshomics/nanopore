# Validation of split basecalling

Rule, scripts and results of the validation of `basecall_split_jobs` (main README, "Split basecalling
across GPU nodes"). Test data: the four smallest pod5 files of run
20261008_1522_P2S-02398-A_PBM28996_eebc9374 (`_77`, `_80`, `_81`, `_82`; 6.47 GB), basecalled with
dorado 2.1.2, `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`, `4mC_5mC@v1`, `6mA@v1` and `--kit-name SQK-NBD114-24`,
as in the run. Every GPU job of the arms ran on standard.q node601–609 (4 × TITAN Xp 12 GB), so the GPU
type cannot explain a difference between arms.

## Rule

`RULES.md` (rule R5: question, arms, five pass criteria, and what is only reported) was written before
any validation job was submitted. Amendment A1 (shipped resources; the incremental re-planning test)
was added before the first submission, and amendment A2 (`h_vmem` 4 × 80G) after that submission
failed and before the second. `results/rules_record.txt` holds each version's sha256 and write time:
R5 at 2026-10-21T04:31:49Z, with A1 at 2026-10-21T04:49:14Z, with A1 and A2 at 2026-10-21T05:07:22Z; the first submission was at 2026-10-21T04:50:34Z and the second at 2026-10-21T05:07:58Z.

## Arms

| arm | pipeline | GPU jobs | hosts | basecalling |
|---|---|---|---|---|
| S1 | `basecall_split_jobs: 1` (the single job) | 1 | node607 | 2 h 13 min (2.9 GB/h) |
| S3 | `basecall_split_jobs: 3`, chunk BAMs kept | 3, then the merge | node606, node602, node604 | 52 min from the first chunk's start to the last chunk's end (chunks 47, 35, 52 min); merge 3 s |
| S1′ | S1 run again (R5 criterion 2, because S3 and S1 differed) | 1 | node603 | 2 h 9 min (3.0 GB/h) |
| P | the records of the four files in the production `basecalls.bam` (`prod_extract.sge`) | – | run 20261012, 2 × RTX A5000 | – |

Criterion 5 used further copies of S3: S3neg (faults b–d, then recovery e), S3size (f), S3inc (the
incremental test g of A1), and S3 itself (re-submission r1).

## Results

| R5 | criterion | result | files in `results/` |
|---|---|---|---|
| 1 | the same read ids, each once | met: 36960 reads in each, none only in one arm, no duplicates | `S3_vs_S1.json` |
| 2 | identical records, or no more differing reads than d′ + 3·√max(1, d′) | 36916 of 36960 records identical; 44 differ only in ML (44) and MM (2): in 42 reads 1–7 per-base modification probabilities differ, by at most 7 of 255. S1′ against S1: 44 reads differ (ML 44, MM 2), by at most 7 of 255; 44 of them are the reads in which S3 differs from S1, and S3 against S1′ differs in 0 (all 36960 records identical). Bound d′ + 3·√max(1, d′) = 63.9: **met** | `S3_vs_S1.json`, `S1b_vs_S1.json`, `inspect_diffs.json`, `p3/` |
| 3 | the merged file is intact and consistent | met: `quickcheck -u` passes; 36960 records = the chunks' 36960 = S1's; the 6 read groups used are all declared (the same 25 `@RG` lines as S1); `@PG` ids unique; 4 `@CO` lines record the plan | `quickcheck.txt`, `S3_merged_header.sam`, `S3_chunk_summary.tsv` |
| 4 | the same reads per barcode after demux | met: `s02_demux` 6 FASTQ files (36960 reads) and `s02b_demux_bam` 6 BAMs (36960 reads), the same read ids per file as from S1 | `demux_S3_vs_S1.json` |
| 5 | failures stop the stage; resume and re-planning | met: (b) a missing chunk, (c) a pod5 file outside the plan and (d) a truncated chunk each stopped the merge with FAIL and no `basecalls.bam`, (d) moving the chunk to `.bad`; (f) a changed file size stopped task 1 before dorado; (e) re-submitting re-basecalled only chunk 01 and restored S1's 36960 reads; (g) after a fifth file was added, re-generating kept chunks 01–03 and added chunk 04, the only chunk then basecalled (36960 reads = S1's, plus the fifth file's 13254, with production's parent reads); (r1) re-submitting the finished S3 skipped every task and the merge | `neg_tests.tsv`, `neg_*.out`, `inc_plan_*.tsv`, `p2/` |
| reported | against the production BAM (other hardware) | R5 asks for the same parent reads: **not met**. 1 parent read is only in S3 (35 bases, mean Q 3.98, also in S1) and 1 only in production (177 bases, Q 7.45, not in S1). S1 differs from production in exactly the same 127 + 101 reads, parent reads and fields, so the hardware causes this, not the split. None of the 36833 shared records is identical across hardware | `S3_vs_P.json`, `S1_vs_P.json`, `inspect_diffs.json` |
| reported | throughput | TITAN Xp: 2.89, 2.93, 2.90 GB/h for the three chunks, 2.9 GB/h for S1; the generator's warning uses 2.8 GB/h | `chunk_0*.done`, `S1_pipeline.log`, `S3_pipeline.log` |

`compare_bams.py` reduces each record to digests of FLAG, SEQ, QUAL and every tag. Its positive control
(`make_control.py`: a copy of S1 without its first read and with one base changed) was reported as
1 read missing and 1 read differing (qs, seq) (`p2/control_S1mod_vs_S1.json`).

## Problems the validation found (fixed before the results above)

1. **`h_vmem` and CUDA start-up.** At amendment A1's 4 × 10G, dorado stopped within seconds in every
   job, the unchanged single job included, with "CUDA out of memory" while most of the card was free.
   Two 3,000-read probes (`probe_dorado.sge`) separated the causes: at 4 × 10G (`ulimit -v` 40 GiB) on node705, dorado exited with code 1 after 9 s ("CUDA out of memory" for 2.67 GiB on a 10.75 GiB card with 6.99 GiB free); at 4 × 80G on node607 (TITAN Xp) it basecalled 3007 reads in 477 s. SGE
   enforces `h_vmem` as the process's address-space limit, and CUDA start-up reserves far more address
   space than dorado then uses. The split jobs therefore request 4 × 80G (`basecall_split_h_vmem`, A2).
2. **`samtools quickcheck` and unaligned BAMs.** Without `-u`, `samtools quickcheck` rejects an intact unaligned BAM (exit 8, "had no targets in header"), so the merge as first written would have rejected every intact chunk (found while diagnosing the first submission). With `-u` the intact BAM passes (exit 0) and a truncated chunk still fails (exit 16, "missing EOF block"). The merge therefore checks chunks with
   `quickcheck -u` (`p3/quickcheck_modes.txt`).
3. **Activating the conda environment in the merge job.** Activating the `nanopore` environment runs
   `R CMD javareconf`, which starts R. The first re-submitted merge job stayed in that step on node514
   (job 6269360: `R --vanilla`, started by `javareconf`, at 05:57 elapsed and 00:06:00 CPU time; `attempt1/bsv_ps514c.log`). In `act_test.sge`, run on node514 with the merge job's
   limits, `mamba activate nanopore` did not finish within 300 s (4 slots, `ulimit -v` 16 GiB, 128 cores); with the BLAS libraries limited to one thread it took 2 s (`attempt1/bsv_act514.log`); on three other nodes
   activation took 0.4–0.6 s (node517, node524, node527) (`probe_activate.sge`, `attempt1/bsv_act_node*.log`). `s01b_merge_bam`
   therefore sets `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS` and `MKL_NUM_THREADS` to 1 before activating.
4. **Holds of the validation jobs.** The generated stages hold on job names, and `qsub -hold_jid` on
   the command line adds to a script's own `-hold_jid` rather than replacing it. The validation jobs
   (named `bsv_…`/`bsv2_…`) were therefore re-pointed with `qalter -hold_jid <job ids>`, so that they
   could neither wait for nor release the jobs of another pipeline run.

## Scripts, in the order they ran

1. `val_setup.sh` – test pod5 folders (symlinks), configs and the generated pipelines S1, S3, S3neg, S3size.
2. `val_submit_main.sh` – S1 and S3 (GPU jobs pinned to node601–609), the S3 merge, demux of both arms,
   arm P (`prod_extract.sge`) and the comparison (`bsv_compare.sge`: `compare_bams.py`, `demux_compare.py`).
3. `val_setup2.sh`, `val_submit_p2.sh`, `val_submit_p2b.sh` – criterion 5: the CPU driver `neg_driver.sge`
   (faults b–d, first half of g), then e, f, r1 and g, and the comparisons (`compare_p2.sge`).
4. `inspect_diffs.sge` (`inspect_diffs.py`) – the parent reads found in only one of S3, S1 and P, and the
   size of the ML/MM differences.
5. `val_submit_s1b.sh` – S1′; then `compare_s1b.sge` – S1′ against S1 and S3, with `inspect_diffs.py` and
   `overlap_diffs.py`. `ovtest.sge` checked `overlap_diffs.py` first, with S1 in place of S1′: it found the
   44 reads of `S3_vs_S1.json` and none between S1 and itself (`p3/selftest_overlap_S3_S1_S1.json`).

`probe_dorado.sge`, `probe_activate.sge` and `act_test.sge` are the diagnostic probes of problems 1 and 3. Every script
writes below `/ebio/abt3_scratch/jmarsh/software/science/split_build/val`; `results/jobs*.tsv` list the
job ids, and `results/qacct_*.txt` the accounting records.
