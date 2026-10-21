# Rule R4: acceptance of the s18e_phava stage (fixed 2026-10-18, before any validation result)

PhaVa (v0.2.3) is run as published with its recommended report filter: an inverton is
**called** when at least 3 reads support the inverted orientation and they are at least 3 % of
the reads PhaVa counts there (forward + reverse). The stage adds: a deterministic read
subsample to `phava_max_depth` (the value fixed by the timing test and written to
`evidence/val/design.json` before any spiked run), a PhaVa-ready GFF, and gene annotation.
R4 tests that this integration measures inverted fractions correctly on real reads.

## Design

- Samples: Cminuta01 (Bakta GFF3) and Tube14_T002_conj (Prokka GFF).
- Reads: all primary and unmapped reads of the sample's s13 mapping, subsampled by read name
  (seeded) to `phava_max_depth`; the same reads in every condition of a sample.
- Control: PhaVa (`variation_wf --reportAll`) on the unmodified reads.
- Spiked invertons: 20 per sample, drawn with seed 1 from the control invertons with 0 reverse
  reads, at least 100 forward reads and a middle segment of 30-2,000 bp, at least 10 kb from
  each other and from every inverton with a reverse read in the control, and overlapping no
  other PhaVa candidate (amendment A1).
- Edit: in a seeded random fraction f of the reads whose primary alignment spans the
  inverton with 50 bp to spare on each side, the read bases aligned to the middle segment
  (between the inverted repeats) are reverse-complemented and their qualities reversed. This
  is the sequence PhaVa uses for the inverted orientation. Fractions f = 0.05, 0.20 and 0.50.
  Truth per inverton: f_real = edited spanning reads / spanning reads.
- Each spiked read set is analysed like the control, in its own PhaVa directory.

## Acceptance (all must hold in both samples)

1. **Detection**: in every spiked condition at least 19 of the 20 spiked invertons are called.
2. **Accuracy**: for the called spiked invertons of a condition, |ProportionReverse - f_real|
   <= 0.05 for at least 90 % of them, and the median absolute difference is <= 0.02.
3. **No spill-over**: in every spiked condition the set of called unspiked invertons equals the
   control's.
4. **Clean control**: no spiked inverton has a reverse read in the control (checks the draw).

If R4 fails, the cause is found and fixed in the integration (read preparation, GFF, parsing,
counting) and the whole validation is re-run; the rule is not changed.

## Amendment A1 (2026-10-18, before any control or spiked PhaVa run)

PhaVa keeps candidates that overlap each other (in the March runs: 9 of 1,390 C. minuta and
241 of 2,012 M. smithii candidates, 51 M. smithii pairs sharing at least 80 % of their
segment). Editing one of them changes the bases of the others, so they are neither spiked nor
unspiked in the sense of criterion 3. The draw therefore takes only candidates that overlap
no other candidate. The acceptance criteria are unchanged. In production every inverton lists
its overlapping candidates (column overlapping_candidates).

## Amendment A2 (2026-10-18, after the two control runs, before any spiked result)

The control runs showed that PhaVa's joint mapping of all candidates loses forward reads where
candidates are close (counted reads: median ~900 beyond 1 kb from the next candidate, 63-102
for overlapping ones), while inverted reads, unique to one reference, are kept. The inverted
fraction is then inflated: Cminuta01 1:246611-246624-246728-246741, PhaVa 23/89 = 25.8 %, the
same reads counted against that inverton's two orientations alone 24/809 = 3.0 % (s18d group
G1: 3.0 %). Production therefore re-counts every candidate with at least one PhaVa inverted
read against its own two orientations with PhaVa's counting rules (phava_tools.py genotype,
same reads), and calls and fractions use these counts (PhaVa's counts are kept as
phava_forward_reads / phava_reverse_reads). Criteria 1-3 are evaluated on the production
values (status and reverse_fraction of invertons.tsv); PhaVa's own ProportionReverse is
reported alongside. The thresholds, fractions, draw and acceptance limits are unchanged.
