#!/usr/bin/env python3
"""Diagnostic: genotype invertons with only their own two orientations as references.

For each inverton (PhaVa id contig:a-b-c-d), every read with an alignment overlapping the inverton
in the s13 BAM is mapped (full read, mappy map-ont) to two templates: the consensus with 1 kb flanks
and the same with the segment b+1..c reverse-complemented. PhaVa's counting rules are applied to the
primary hit: MAPQ > 2, at most 15 % mismatches+indels over the aligned part (NM / block length), the
hit spans the inverton with more than 20 bp on each side. A read hashing outside the subsample
fraction (--fraction, seed phava1) is skipped so the same reads as the PhaVa run are used.
argv: --bam --consensus --invertons (ids, one per line) --out"""
import argparse
import hashlib
import math
import os
import tempfile

import mappy
import pysam

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(s):
    return s.translate(COMP)[::-1]


def hash_unit(seed, rid):
    return int.from_bytes(hashlib.blake2b(f'{seed}:{rid}'.encode(), digest_size=8).digest(), 'big') / 18446744073709551616.0


def wilson(k, n, z=1.96):
    if not n:
        return float('nan'), float('nan')
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, (c - h) / d), min(1.0, (c + h) / d)


ap = argparse.ArgumentParser()
ap.add_argument('--bam', required=True)
ap.add_argument('--consensus', required=True)
ap.add_argument('--invertons', required=True)
ap.add_argument('--fraction', type=float, default=1.0)
ap.add_argument('--flank', type=int, default=1000)
ap.add_argument('--out', required=True)
a = ap.parse_args()
cons = {r.name: r.sequence.upper() for r in pysam.FastxFile(a.consensus)}
bam = pysam.AlignmentFile(a.bam)
out = open(a.out, 'w')
out.write('inverton\treads_tested\tforward\treverse\tunassigned\tfraction\tci_low\tci_high\n')
for inv in [l.strip() for l in open(a.invertons) if l.strip()]:
    ctg, co = inv.rsplit(':', 1)
    A, B, C, D = map(int, co.split('-'))
    s = cons[ctg]
    lo, hi = max(0, A - a.flank), min(len(s), D + a.flank)
    F = s[lo:hi]
    R = s[lo:B] + rc(s[B:C]) + s[C:hi]
    i0, i1 = A - lo, D - lo                     # inverton within the templates
    with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as fh:
        fh.write(f'>f\n{F}\n>r\n{R}\n')
        fa = fh.name
    al = mappy.Aligner(fa, preset='map-ont')
    os.unlink(fa)
    seqs = {}
    for x in bam.fetch(ctg, A, D):
        if x.is_secondary or x.query_name in seqs:
            continue
        if a.fraction < 1 and hash_unit('phava1', x.query_name) >= a.fraction:
            continue
        if x.is_supplementary:
            seqs.setdefault(x.query_name, None)
            continue
        seqs[x.query_name] = x.query_sequence
    missing = [n for n, q in seqs.items() if q is None]
    if missing:   # reads seen only as supplementary here: take the primary record's sequence
        want = set(missing)
        for x in bam.fetch(ctg, max(0, A - 60000), min(len(s), D + 60000)):
            if x.query_name in want and not x.is_secondary and not x.is_supplementary:
                seqs[x.query_name] = x.query_sequence
                want.discard(x.query_name)
    nf = nr = nu = 0
    for name, q in seqs.items():
        if q is None:
            nu += 1
            continue
        hits = [h for h in al.map(q) if h.is_primary]
        if not hits:
            nu += 1
            continue
        h = hits[0]
        span = h.r_st < i0 - 20 and h.r_en > i1 + 20
        if h.mapq > 2 and span and h.NM <= 0.15 * h.blen:
            if h.ctg == 'f':
                nf += 1
            else:
                nr += 1
        else:
            nu += 1
    n = nf + nr
    l, u = wilson(nr, n)
    out.write(f'{inv}\t{len(seqs)}\t{nf}\t{nr}\t{nu}\t{nr / n if n else float("nan"):.4f}\t{l:.4f}\t{u:.4f}\n')
    out.flush()
out.close()
