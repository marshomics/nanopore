# Pre-specified decision rules (written 2026-10-17, before any validation result)

These rules decide stage defaults. They are fixed before any validation result
is seen. A held-out result never feeds back into a choice.

## R1 — consensus polishing (s12b_polish)

For each sample s and consensus state x (draft = Autocycler consensus;
polished@d = `dorado polish --bacteria` on a deterministic d-fold read subsample,
d in {50, 100, 200, 400}):

- Dref(s,x): SNPs + indels of at most 50 bp between x and the independent
  reference (C. minuta: NCBI NZ_CP029256; M. smithii: the user's T002_i2),
  counted from a minimap2 asm5 alignment (cs tag). Larger events are counted
  separately and do not enter Dref.
- Drep(p,x): the same count between the two replicates of pair p
  (Cminuta01 vs Cminuta02; Tube14 vs Tube15), both in state x.
- H(s,x): PASS calls from Clair3 (haploid_precise, model r1041_e82_400bps_sup_v520)
  when a 100-fold subsample of hold-out reads is mapped to x. Hold-out reads are
  disjoint from every polishing subsample, split by read-name hash.
- E(s,d): of 300 injected SNPs and 300 injected homopolymer +-1 indels in the
  polished consensus, the fraction polishing at depth d restores, and the number
  of new differences it introduces. Reported, not used for the decision.

Default depth d* = the smallest tested d with
  sum_s Dref(s,d) <= min_d' sum_s Dref(s,d') + 1  and  sum_s H(s,d) <= min_d' sum_s H(s,d') + 2.

polish_enabled defaults to true if and only if, at d*:
  sum_s Dref(polished) <= sum_s Dref(draft),
  sum_p Drep(polished) <= sum_p Drep(draft),
  sum_s H(polished)    <= sum_s H(draft).
Otherwise it defaults to false. The sums run over all 4 samples and both pairs.
If the joint rule holds but fails for the two M. smithii samples alone, the README
says so explicitly.

## R2 — small variants between samples and against references (s18c_snv)

Truth comes from a mutated reference. Each sample's polished consensus gets 1,500
SNPs and 300 indels (1-20 bp; half change a homopolymer of at least 3 bp by +-1-2 bp),
at uniform random positions at least 50 bp apart and at least 500 bp from contig
ends; 10 % of them sit in repeat-masked regions. The sample's real reads, as
deterministic subsamples at 10, 25, 50 and 100-fold, are mapped with minimap2 -x lr:hq.

Candidate configurations:
- C3P: Clair3 2.x --haploid_precise
- C3AF: Clair3 diploid output; a het becomes hom-ALT if AF >= 0.65, otherwise it is dropped
- DP: dorado polish --bacteria --vcf
- INT: C3AF intersect DP (normalised alleles)
- UNI: C3AF union DP

Scoring uses vcfdist on PASS calls, giving F1 for SNPs and indels combined.
The tuning set is Cminuta01 and Cminuta02; the held-out set is Tube14 and Tube15.
The chosen configuration has the highest pooled tuning F1 at 50-fold. A tie
(|dF1| < 0.0005) goes to higher precision, then to a single caller over a
combination. The snv_depth default is the smallest depth in {10, 25, 50, 100} whose
pooled tuning F1 for the chosen configuration is at least F1(100-fold) - 0.001.
Held-out metrics are reported for every configuration and depth, and never used
for selection.

## R3 — minority variants within a sample (s18d_minor)

Engine defaults are fixed before validation: FDR q <= 0.01 (Benjamini-Hochberg),
at least 3 ALT reads on each strand, strand depth at least 20, base quality at
least 10, MAPQ at least 20, and a repeat mask (self-alignment of at least 300 bp
at 85 % identity or more).

Spike-in design: real reads of Cminuta01 and Tube14 at 100, 500 and 2,000-fold
(deterministic read subsample). A minor haplotype is a random set of reads with
fraction f in {0.25, 0.5, 1, 2, 5, 10, 20, 40} %. At 300 sites (200 isolated; 20
clusters of 5 sites 1 kb apart) plus 30 sites inside masked repeats, the minor reads'
base is replaced by a random non-reference base. The truth fraction is the realised
share of minor reads among the reads counted at that site.

Acceptance, per condition (sample x depth x fraction):
- False-discovery proportion. A false call is a call at a non-spiked, non-masked
  position that is not also called in the unmodified control of the same sample
  and depth. With 20 or more calls, false calls / all calls must be <= 0.05;
  with fewer than 20, false calls must be <= 1.
  Calls present in controls are reported separately as background and characterised.
- Among called spiked sites with truth fraction >= 1 %, the reported 95 % interval
  contains the truth for >= 90 %.
- No spiked site inside a masked repeat is called; these sites are reported as not assessed.

Sensitivity, lowest detected fraction and fraction error are measured and
reported, not thresholded.
