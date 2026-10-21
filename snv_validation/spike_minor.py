#!/usr/bin/env python3
"""spike_minor.py - spike a minor haplotype into real reads (validation of s18d_minor only).

  sites   choose spike-in sites from the engine's repeat mask: isolated sites, clusters of
          sites 1 kb apart (testable for linkage) and sites inside masked repeats
  spike   copy a BAM; reads whose name hashes below --fraction form the minor haplotype and
          get the ALT base at every site they cover. Writes the truth table: the realised
          share of minor reads among reads the engine counts at each site (MAPQ >= 20,
          base quality >= 10, primary or supplementary, not duplicate/QC-fail).

Only aligned bases (CIGAR M/=/X) are edited; read length, CIGAR and qualities are unchanged,
so the reads keep their real error profile everywhere else, including at the site in every
read of the major haplotype. MD/NM tags are removed because they would be stale.
"""
import argparse
import bisect
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import minor_variants as mv      # noqa: E402  (same read filter and hashing as the engine)


def cmd_sites(a):
    import pysam
    rng = random.Random(a.seed)
    recs = mv.read_fasta(a.ref)
    reps = mv.load_bed(a.repeats)
    seqs = {n: s for n, s, _ in recs}
    order = [n for n, _, _ in recs]
    total = sum(len(s) for s in seqs.values())
    bam = pysam.AlignmentFile(a.bam)
    taken = {n: [] for n in order}

    def masked(n, i):
        return any(s <= i < e for s, e in reps.get(n, []))

    def free(n, s, e, gap):
        return all(e + gap <= s0 or s - gap >= e0 for s0, e0 in taken[n])

    def rand_site(want_masked, gap, span=1):
        for _ in range(500000):
            x = rng.randrange(total)
            for n in order:
                if x < len(seqs[n]):
                    break
                x -= len(seqs[n])
            if x < a.end_mask + 50 or x + span >= len(seqs[n]) - a.end_mask - 50:
                continue
            if any(seqs[n][x + k * a.cluster_spacing] not in 'ACGT' for k in range(span)):
                continue
            if want_masked != masked(n, x):
                continue
            if not want_masked and any(masked(n, x + k * a.cluster_spacing) for k in range(span)):
                continue
            if not free(n, x, x + (span - 1) * a.cluster_spacing + 1, gap):
                continue
            return n, x
        raise SystemExit('could not place spike-in sites')

    out = []
    for _ in range(a.n_clusters):
        n, x = rand_site(False, a.min_gap, a.cluster_size)
        taken[n].append((x, x + (a.cluster_size - 1) * a.cluster_spacing + 1))
        for k in range(a.cluster_size):
            p = x + k * a.cluster_spacing
            out.append((n, p, 'cluster'))
    for _ in range(a.n_isolated):
        n, x = rand_site(False, a.min_gap)
        taken[n].append((x, x + 1))
        out.append((n, x, 'isolated'))
    for _ in range(a.n_masked):
        n, x = rand_site(True, 500)
        taken[n].append((x, x + 1))
        out.append((n, x, 'masked'))
    with open(a.out, 'w') as fh:
        fh.write('contig\tpos\tref\talt\tset\n')
        for n, p, st in sorted(out, key=lambda t: (order.index(t[0]), t[1])):
            ref = seqs[n][p]
            alt = rng.choice([b for b in 'ACGT' if b != ref])
            fh.write(f'{n}\t{p + 1}\t{ref}\t{alt}\t{st}\n')
    print(f'{len(out)} sites written to {a.out}')


def aligned_qpos(read, refpos):
    """Query index aligned (M/=/X) to reference position refpos, or None."""
    rpos, qpos = read.reference_start, 0
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):
            if rpos <= refpos < rpos + ln:
                return qpos + (refpos - rpos)
            rpos += ln
            qpos += ln
        elif op in (1, 4):
            qpos += ln
        elif op in (2, 3):
            if rpos <= refpos < rpos + ln:
                return None
            rpos += ln
        if rpos > refpos:
            return None
    return None


def cmd_spike(a):
    import pysam
    sites = [l.rstrip('\n').split('\t') for l in open(a.sites)][1:]
    by = {}
    for ctg, pos, ref, alt, st in sites:
        by.setdefault(ctg, []).append((int(pos) - 1, ref, alt, st))
    for v in by.values():
        v.sort()
    keys = {c: [s[0] for s in v] for c, v in by.items()}
    tally = {(c, s[0]): [0, 0] for c, v in by.items() for s in v}      # [minor counted, all counted]
    flt = mv.ReadFilter(20, 1.0, 'none')
    inp = pysam.AlignmentFile(a.bam)
    out = pysam.AlignmentFile(a.out_bam, 'wb', template=inp)
    n_reads = n_minor = n_edits = 0
    for r in inp.fetch(until_eof=True):
        n_reads += 1
        if r.is_unmapped or r.reference_name not in by or r.cigartuples is None or r.query_sequence is None:
            out.write(r)
            continue
        ks = keys[r.reference_name]
        lo = bisect.bisect_left(ks, r.reference_start)
        hi = bisect.bisect_left(ks, r.reference_end)
        if lo == hi:
            out.write(r)
            continue
        minor = mv.name_unit(r.query_name, a.seed + ':minor') < a.fraction
        counted = flt.ok(r)
        quals = r.query_qualities
        seq = None
        for p, ref, alt, st in by[r.reference_name][lo:hi]:
            q = aligned_qpos(r, p)
            if q is None:
                continue
            if counted and quals is not None and quals[q] >= 10:
                t = tally[(r.reference_name, p)]
                t[1] += 1
                t[0] += minor
            if minor:
                if seq is None:
                    seq = list(r.query_sequence)
                seq[q] = alt
                n_edits += 1
        if seq is not None:
            r.query_sequence = ''.join(seq)
            r.query_qualities = quals
            for tag in ('MD', 'NM'):
                if r.has_tag(tag):
                    r.set_tag(tag, None)
            n_minor += 1
        out.write(r)
    out.close()
    pysam.index(a.out_bam)
    with open(a.out_truth, 'w') as fh:
        fh.write('contig\tpos\tref\talt\tset\ttruth_fraction\ttruth_minor_reads\ttruth_depth\n')
        for ctg, pos, ref, alt, st in sites:
            m, n = tally[(ctg, int(pos) - 1)]
            fh.write(f'{ctg}\t{pos}\t{ref}\t{alt}\t{st}\t{m / n if n else 0:.6f}\t{m}\t{n}\n')
    print(f'{n_reads} reads, {n_minor} edited minor reads, {n_edits} base edits, fraction {a.fraction}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('sites')
    p.add_argument('--ref', required=True); p.add_argument('--repeats', required=True)
    p.add_argument('--bam', required=True); p.add_argument('--out', required=True)
    p.add_argument('--n-isolated', type=int, default=200); p.add_argument('--n-clusters', type=int, default=20)
    p.add_argument('--cluster-size', type=int, default=5); p.add_argument('--cluster-spacing', type=int, default=1000)
    p.add_argument('--n-masked', type=int, default=30); p.add_argument('--min-gap', type=int, default=2000)
    p.add_argument('--end-mask', type=int, default=500); p.add_argument('--seed', type=int, default=17)
    p = sp.add_parser('spike')
    p.add_argument('--bam', required=True); p.add_argument('--sites', required=True)
    p.add_argument('--fraction', type=float, required=True); p.add_argument('--seed', default='spike1')
    p.add_argument('--out-bam', required=True); p.add_argument('--out-truth', required=True)
    a = ap.parse_args()
    {'sites': cmd_sites, 'spike': cmd_spike}[a.cmd](a)


if __name__ == '__main__':
    main()
