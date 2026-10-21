# Rule R5: split basecalling (s01_basecall array + s01b_merge_bam)

Written 2026-10-21, before any validation job was submitted. The write time and sha256 of this file
are recorded in `results/rules_record.txt`.

## Question

Does basecalling the pod5 files in several GPU jobs and concatenating the chunk BAMs give the same
`basecalled/basecalls.bam` content, and the same demultiplexed reads, as one dorado job over the
whole folder?

## Test data

The four smallest pod5 files of run `20261008_1522_P2S-02398-A_PBM28996_eebc9374`
(`PBM28996_eebc9374_9cc5eed7_77/80/81/82.pod5`, 6.47 GB), symlinked into one test folder. The
basecalling settings are those of the run: `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`, the `4mC_5mC` and
`6mA` models, `--kit-name SQK-NBD114-24`, dorado 2.1.2.

## Arms

- **S1**: the pipeline generated with `basecall_split_jobs: 1` (the existing single-job `s01_basecall`).
- **S3**: the pipeline generated with `basecall_split_jobs: 3` (three array tasks, then
  `s01b_merge_bam`).
- Both arms' GPU jobs are restricted to the same node class (standard.q node601–node609, 12 GB cards),
  so that a difference cannot come from GPU type.
- **P** (reference, not a pass criterion except read coverage): the records of the production
  `basecalls.bam` of this run whose `fn` tag names one of the four files (one job, other hardware).

## Pass criteria (all required)

1. **Reads.** S3 and S1 hold the same read ids, each exactly once.
2. **Records.** For every read, FLAG, SEQ, QUAL and every tag (including MM, ML, RG and the barcode
   tags) are identical between S3 and S1. If any read differs, a second single-job arm S1' is run
   and S3 passes only if the number of reads differing between S3 and S1 is at most
   d' + 3·sqrt(max(1, d')), where d' is the number differing between S1' and S1.
3. **Merged file.** `samtools quickcheck` passes; the record count equals the sum of the chunk counts
   and the S1 count; every RG id used by a record is declared in the header; @PG ids are unique.
4. **Demux.** `s02_demux` (FASTQ) and `s02b_demux_bam` (BAM) run on S3 give, per barcode, the same
   read ids as on S1.
5. **Failure handling.** Each of these stops the stage with a FAIL event in `pipeline.log` and no
   `basecalls.bam`: a missing chunk, a corrupt chunk, a pod5 file in the folder that the manifest
   lacks; a pod5 file whose size differs from the manifest stops its task before dorado starts, with
   no chunk BAM. Re-submitting after a completed run skips every task and the merge. Re-submitting
   after a lost chunk re-basecalls only that chunk, and the merged read set again equals S1's.

## Reported, not pass criteria

- S3 against P: parent-read coverage (required: the same parent reads) and the fraction of reads
  whose records are identical across GPU hardware.
- Throughput of each GPU job (GB of pod5 per hour, by node and card), used to choose the shipped
  number of jobs and stated in the README.

## Amendment A1 (2026-10-21, before any validation job was submitted)

- The split jobs request the shipped resources, 4 slots × 10G (40 GB of virtual memory) instead of
  8 slots × 8G: the 12 GB-card nodes rarely have more than a few free slots, and dorado has used
  4.1–4.6 cores and at most 26 GB of virtual memory in full runs here. Arm S1 uses the same resources.
- Re-generating into the same output folder now keeps the finished chunks when pod5 files were added
  and plans only the new files into extra chunks. Added to criterion 5: after a fifth pod5 file
  (`PBM28996_eebc9374_9cc5eed7_79.pod5`) is added to the folder of a finished split run and the
  pipeline is re-generated, only the new chunk is basecalled, and the merged parent reads equal those
  of S1 plus the production (P) parent reads of the fifth file.

## Amendment A2 (2026-10-21, after the first submission failed and before the second)

- At A1's 4 slots × 10G, dorado stopped within seconds of starting, in arm S1 (the unchanged
  single-job script) and in every S3 task, with "CUDA out of memory" while most of the card was
  free. Two 3,000-read probes separated the causes: at 4 × 80G a TITAN Xp node (node601–609)
  basecalled normally, and at 4 × 10G an 11 GB-card node (node705) failed the same way. The limit
  is SGE's h_vmem (enforced as the process's address-space limit), not the card type: dorado's CUDA
  start-up needs far more address space than the ≤26 GB it then uses. The split jobs and arm S1
  therefore request 4 slots × 80G (320 GB), the per-slot value the single job already used. Nothing
  else in the rule changes; the failed first submission is not part of the results.
