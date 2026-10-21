#!/usr/bin/env python3
"""spike_inv.py - read-level inversion spike-ins for the s18e_phava validation (rule R4).

  select  draw the spiked invertons from a control inverton table (phava_tools annotate):
          0 reverse reads, >= --min-forward forward reads, segment --min-seg..--max-seg bp,
          no overlapping candidate (amendment A1), at least --min-gap bp from each other and
          from every inverton with a reverse read in the control; seeded shuffle, greedy
  edit    in a seeded fraction of the reads whose primary alignment (s13 BAM) spans a spiked
          inverton with --spare bp on each side, reverse-complement the read bases aligned to
          the segment between the inverted repeats and reverse their qualities; one FASTQ per
          fraction plus the truth table (spanning and edited reads per inverton)

Reads are edited in their FASTQ orientation: for a reverse-strand alignment the segment's
query interval is mirrored onto the original read. The edited fractions are nested (the
reads edited at 5 % are among those edited at 20 %).
"""
import argparse
import bisect
import json
import random
import sys

import pysam

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(s):
    return s.translate(COMP)[::-1]


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        return [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in fh if l.strip()]


def gap(a, b, clen):
    """Bases between spans a=(s, e) and b on a circular contig of length clen (0 if they overlap)."""
    lin = max(b[0] - a[1], a[0] - b[1], 0)
    wrap = clen - max(a[1], b[1]) + min(a[0], b[0])
    return min(lin, max(wrap, 0)) if clen else lin


def cmd_select(a):
    rows = read_tsv(a.control)
    clen = {}
    for name, seq in ((r.name, r.sequence) for r in pysam.FastxFile(a.genome)):
        clen[name] = len(seq)
    rev = [r for r in rows if int(r['reverse_reads']) > 0]
    elig = []
    for r in rows:
        if (int(r['reverse_reads']) == 0 and int(r['forward_reads']) >= a.min_forward
                and a.min_seg <= int(r['segment_bp']) <= a.max_seg and not r['overlapping_candidates']):
            sp = (int(r['start']), int(r['end']))
            if all(x['contig'] != r['contig'] or gap(sp, (int(x['start']), int(x['end'])), clen[r['contig']]) >= a.min_gap
                   for x in rev):
                elig.append(r)
    elig.sort(key=lambda r: (r['contig'], int(r['start'])))
    rng = random.Random(a.seed)
    rng.shuffle(elig)
    chosen = []
    for r in elig:
        sp = (int(r['start']), int(r['end']))
        if all(x['contig'] != r['contig'] or gap(sp, (int(x['start']), int(x['end'])), clen[r['contig']]) >= a.min_gap
               for x in chosen):
            chosen.append(r)
        if len(chosen) == a.n:
            break
    if len(chosen) < a.n:
        raise SystemExit(f'only {len(chosen)} eligible invertons (of {len(elig)} passing the filters)')
    chosen.sort(key=lambda r: (r['contig'], int(r['start'])))
    cols = ['inverton', 'contig', 'start', 'end', 'segment', 'segment_bp', 'ir_bp', 'forward_reads', 'class']
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in chosen:
            fh.write('\t'.join(r[c] for c in cols) + '\n')
    print(json.dumps({'candidates': len(rows), 'with_reverse_reads': len(rev), 'eligible': len(elig),
                      'chosen': len(chosen)}))


def segment_query_interval(aln, s0, e0):
    """Query interval (BAM orientation, inclusive) aligned to reference s0..e0 (0-based, inclusive)."""
    pairs = aln.get_aligned_pairs(matches_only=True)
    rpos = [p[1] for p in pairs]
    i = bisect.bisect_left(rpos, s0)
    j = bisect.bisect_right(rpos, e0) - 1
    if i >= len(pairs) or j < 0 or i > j:
        return None
    return pairs[i][0], pairs[j][0]


def cmd_edit(a):
    names = set()
    with open(a.fastq) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                names.add(line[1:].split(None, 1)[0])
    spiked = read_tsv(a.spiked)
    fracs = [float(x) for x in a.fractions.split(',')]
    bam = pysam.AlignmentFile(a.bam)
    edits = {f: {} for f in fracs}          # fraction -> read -> [(qa, qb)] in FASTQ orientation
    truth = []
    for r in spiked:
        start, end = int(r['start']), int(r['end'])
        S, E = (int(x) for x in r['segment'].split('-'))
        span = []
        for aln in bam.fetch(r['contig'], S - 1, E):
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped or aln.query_name not in names:
                continue
            if aln.reference_start > start - 1 - a.spare or aln.reference_end < end + a.spare:
                continue
            if any(op == 5 for op, _ in aln.cigartuples):
                continue
            iv = segment_query_interval(aln, S - 1, E - 1)
            if iv is None:
                continue
            L = aln.query_length
            qa, qb = iv if not aln.is_reverse else (L - 1 - iv[1], L - 1 - iv[0])
            seg = aln.query_sequence[iv[0]:iv[1] + 1]
            span.append((aln.query_name, qa, qb, seg if not aln.is_reverse else rc(seg), L))
        span.sort()
        order = list(span)
        random.Random(f'{a.seed}:{r["inverton"]}').shuffle(order)
        for f in fracs:
            k = round(f * len(order))
            for name, qa, qb, seg, L in order[:k]:
                edits[f].setdefault(name, []).append((qa, qb, seg, L))
            truth.append({'fraction': f, 'inverton': r['inverton'], 'contig': r['contig'], 'start': start, 'end': end,
                          'segment': r['segment'], 'spanning_reads': len(order), 'edited_reads': k,
                          'f_real': f'{k / len(order):.6f}' if order else 'NA'})
    for f in fracs:
        out = f'{a.out_prefix}.f{f:g}.fastq'
        n_edit = 0
        with open(a.fastq) as fh, open(out, 'w') as fo:
            while True:
                h = fh.readline()
                if not h:
                    break
                s, p, q = fh.readline(), fh.readline(), fh.readline()
                name = h[1:].split(None, 1)[0]
                ivs = edits[f].get(name)
                if ivs:
                    seq, qual = s.rstrip('\n'), q.rstrip('\n')
                    for qa, qb, seg, L in ivs:
                        if len(seq) != L or seq[qa:qb + 1] != seg:
                            raise SystemExit(f'{name}: FASTQ read does not match its BAM record at {qa}-{qb}')
                        seq = seq[:qa] + rc(seq[qa:qb + 1]) + seq[qb + 1:]
                        qual = qual[:qa] + qual[qa:qb + 1][::-1] + qual[qb + 1:]
                    s, q = seq + '\n', qual + '\n'
                    n_edit += 1
                fo.write(h + s + p + (q if q.endswith('\n') else q + '\n'))
        need = len(edits[f])
        if n_edit != need:
            raise SystemExit(f'fraction {f}: edited {n_edit} reads, expected {need}')
        print(json.dumps({'fraction': f, 'reads_edited': n_edit, 'fastq': out}))
    cols = ['fraction', 'inverton', 'contig', 'start', 'end', 'segment', 'spanning_reads', 'edited_reads', 'f_real']
    with open(a.truth, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for t in truth:
            fh.write('\t'.join(str(t[c]) for c in cols) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('select')
    p.add_argument('--control', required=True)
    p.add_argument('--genome', required=True)
    p.add_argument('--n', type=int, default=20)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--min-forward', type=int, default=100)
    p.add_argument('--min-seg', type=int, default=30)
    p.add_argument('--max-seg', type=int, default=2000)
    p.add_argument('--min-gap', type=int, default=10000)
    p.add_argument('--out', required=True)
    p.set_defaults(func=cmd_select)
    p = sp.add_parser('edit')
    p.add_argument('--bam', required=True)
    p.add_argument('--fastq', required=True)
    p.add_argument('--spiked', required=True)
    p.add_argument('--fractions', default='0.05,0.2,0.5')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--spare', type=int, default=50)
    p.add_argument('--out-prefix', required=True)
    p.add_argument('--truth', required=True)
    p.set_defaults(func=cmd_edit)
    a = ap.parse_args()
    a.func(a)


if __name__ == '__main__':
    main()
