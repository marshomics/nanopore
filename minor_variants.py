#!/usr/bin/env python3
"""minor_variants.py - minority single-nucleotide variants within a sample (s18d_minor).

Model
  * Every read base with base quality >= --min-baseq on a read with MAPQ >= --min-mapq
    (primary or supplementary, not secondary/duplicate/QC-fail) is counted per strand.
  * The error null for each (strand, reference context, ALT base) class is a beta-binomial
    fitted by moments to the same sample's own positions. Contexts are 7-mers centred on
    the position; sparse classes fall back to 5-mers, 3-mers and finally the base itself.
    Positions that look like real variants are removed from the fit iteratively.
  * A site is tested per strand (P[X >= x]); its p-value is the larger of the two
    (it must be significant on both strands). Benjamini-Hochberg runs over all assessed
    positions x 3 ALT bases. A call needs q <= --fdr and >= --min-alt-per-strand ALT
    reads on each strand.
  * Repeats (self-alignment >= 300 bp at >= 85 % identity, as in s18b_sv) and contig
    ends are masked: those positions are not assessed.
  * Fractions are error-corrected, (x/n - e) / (1 - e), with e the null mean for the site,
    and carry Wilson 95 % intervals transformed the same way.
  * Linkage: candidate sites (called, or p <= --candidate-p with >= 2 ALT reads per strand)
    within --link-max-dist are tested for ALT-ALT co-occurrence on the same reads (one-sided
    Fisher test, odds ratio >= --link-min-or, p <= --link-p). Linked sites form groups. A group
    is reported only if it contains at least one called site, and its sub-threshold members
    are then reported as 'linked'. The group fraction counts reads covering >= 2 group
    sites, assigned minor/major by a 2/3 majority of the sites they cover.
  * Known versions (--version FASTA): SNP differences between each supplied genome version
    and the consensus are genotyped at read level and summarised per version.
"""
import argparse
import hashlib
import json
import math
import os
import re
import sys
import time

import numpy as np

VERSION = '1.1'
BASES = 'ACGT'
BIDX = {b: i for i, b in enumerate(BASES)}
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def rc(s):
    return s.translate(COMP)[::-1]


# ------------------------------------------------------------------------- io
def read_fasta(path):
    recs, name, buf, hdr = [], None, [], ''
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if name is not None:
                    recs.append((name, ''.join(buf).upper(), hdr))
                hdr = line[1:]
                name, buf = hdr.split()[0], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf).upper(), hdr))
    if not recs:
        raise SystemExit(f'no sequences in {path}')
    return recs


def is_circular(header, length):
    h = header.lower()
    if re.search(r'circular=(true|yes|y|1)\b', h):
        return True
    if re.search(r'circular=(false|no|n|0)\b', h):
        return False
    return length > 1_000_000


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def name_unit(name, seed):
    d = hashlib.blake2b(f'{seed}:{name}'.encode(), digest_size=8).digest()
    return int.from_bytes(d, 'big') / 18446744073709551616.0


# -------------------------------------------------------------------- repeats
def merge_intervals(ivs):
    out = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def self_repeats(seq, chunk=1000, step=500, min_len=300, min_ident=0.85):
    import mappy
    if len(seq) < chunk:
        return []
    al = mappy.Aligner(seq=seq, preset='asm20', best_n=50)
    ivs = []
    last = len(seq) - chunk
    starts = list(range(0, last + 1, step))
    if starts[-1] != last:
        starts.append(last)
    for st in starts:
        for h in al.map(seq[st:st + chunk]):
            if h.strand == 1 and abs(h.r_st - (st + h.q_st)) <= 50:
                continue
            if h.blen < min_len or h.mlen < min_ident * h.blen:
                continue
            ivs.append((st + h.q_st, st + h.q_en))
    return merge_intervals(ivs)


def repeat_intervals(recs):
    out = {}
    for name, seq, hdr in recs:
        if is_circular(hdr, len(seq)) and len(seq) > 4000:
            raw = self_repeats(seq + seq[:2000])
            out[name] = merge_intervals([(s, min(e, len(seq))) for s, e in raw if s < len(seq)])
        else:
            out[name] = self_repeats(seq)
    return out


def load_bed(path):
    out = {}
    for line in open(path):
        f = line.split()
        if len(f) >= 3 and not line.startswith(('#', 'track')):
            out.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return {k: merge_intervals(v) for k, v in out.items()}


# ----------------------------------------------------------------- annotation
CODON = {}
_b, _aa = 'TCAG', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b):
            CODON[x + y + z] = _aa[16 * i + 4 * j + k]
STARTS = {'ATG', 'GTG', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA'}


def gff_name_map(ref_recs, fna_path):
    """Annotation contig name -> consensus contig name, by identical sequence (Bakta renames contigs)."""
    if not fna_path or not os.path.exists(fna_path):
        return {}
    by_seq = {hashlib.sha1(seq.encode()).hexdigest(): name for name, seq, _ in ref_recs}
    out = {}
    for n, sq, _ in read_fasta(fna_path):
        h = hashlib.sha1(sq.encode()).hexdigest()
        if h in by_seq:
            out[n] = by_seq[h]
    return out


def load_gff(path, name_map=None):
    feats = {}
    name_map = name_map or {}
    if not path or not os.path.exists(path):
        return feats
    dec = lambda s: re.sub(r'%([0-9A-Fa-f]{2})', lambda m: chr(int(m.group(1), 16)), s)
    with open(path) as fh:
        for line in fh:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] in ('gene', 'region', 'source', 'repeat_region'):
                continue
            attr = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
            feats.setdefault(name_map.get(f[0], f[0]), []).append((int(f[3]), int(f[4]), f[6], f[2],
                                               dec(attr.get('locus_tag', attr.get('ID', ''))),
                                               dec(attr.get('gene', '')), dec(attr.get('product', ''))))
    for v in feats.values():
        v.sort()
    return feats


def snv_effect(seq, feats, pos1, ref, alt):
    hits = [f for f in feats if f[0] <= pos1 <= f[1]]
    if not hits:
        near = min(feats, key=lambda f: min(abs(f[0] - pos1), abs(f[1] - pos1)), default=None)
        if near is None:
            return '', '', '', 'intergenic', ''
        d = min(abs(near[0] - pos1), abs(near[1] - pos1))
        return near[4], near[5], near[6], f'intergenic ({d} bp from {near[4] or near[5]})', ''
    start, end, strand, ftype, locus, gene, product = hits[0]
    if ftype != 'CDS':
        return locus, gene, product, ftype, ''
    cds = seq[start - 1:end]
    new = seq[:pos1 - 1] + alt + seq[pos1:]
    ncds = new[start - 1:end]
    if strand == '-':
        cds, ncds = rc(cds), rc(ncds)
    tr = lambda s: ''.join(CODON.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))
    aa0, aa1 = tr(cds), tr(ncds)
    if cds[:3] in STARTS:
        aa0 = 'M' + aa0[1:]
    if ncds[:3] in STARTS:
        aa1 = 'M' + aa1[1:]
    ch = [(i, x, y) for i, (x, y) in enumerate(zip(aa0, aa1)) if x != y]
    if not ch:
        if cds[:3] != ncds[:3]:
            return locus, gene, product, 'start_codon_changed', f'{cds[:3]}>{ncds[:3]}'
        return locus, gene, product, 'synonymous', ''
    i, x, y = ch[0]
    hg = f'p.{x}{i + 1}{y}'
    if i == 0 and cds[:3] in STARTS and ncds[:3] not in STARTS:
        return locus, gene, product, 'start_lost', hg
    if y == '*':
        return locus, gene, product, 'stop_gained', hg
    if x == '*':
        return locus, gene, product, 'stop_lost', hg
    return locus, gene, product, 'missense', hg


def group_genes(lines):
    """genes column of groups.tsv, from the annotated site rows: per group, each gene in order
    of its first site, with the number of sites per effect class (an intergenic site counts for
    its nearest gene)."""
    per = {}
    for r in lines:
        gid = str(r[18])
        if gid == '-':
            continue
        lt, gene, prod, eff = (str(x).replace('\t', ' ') for x in r[19:23])
        cls = 'intergenic' if eff.startswith('intergenic') else (eff or 'unannotated')
        label = f"{lt}{'(' + gene + ')' if gene else ''}:{prod}" if lt else 'unannotated'
        cnt = per.setdefault(gid, {}).setdefault(label, {})
        cnt[cls] = cnt.get(cls, 0) + 1
    return {gid: ' | '.join(f"{lab} [{', '.join(f'{c} x{n}' for c, n in cls.items())}]" for lab, cls in d.items())
            for gid, d in per.items()}


# ------------------------------------------------------------------- counting
class ReadFilter:
    def __init__(self, min_mapq, keep_frac, seed):
        self.min_mapq, self.keep_frac, self.seed = min_mapq, keep_frac, seed

    def ok(self, r):
        if r.is_unmapped or r.is_secondary or r.is_qcfail or r.is_duplicate:
            return False
        if r.mapping_quality < self.min_mapq:
            return False
        if self.keep_frac < 1.0 and name_unit(r.query_name, self.seed) >= self.keep_frac:
            return False
        return True


def mean_depth(bam, recs, flt):
    tot = 0
    for name, seq, _ in recs:
        for r in bam.fetch(name):
            if flt.ok(r):
                tot += r.reference_length or 0
    return tot / sum(len(s) for _, s, _ in recs)


LUT = np.full(256, 4, dtype=np.int64)
for _i, _b in enumerate('ACGT'):
    LUT[ord(_b)] = _i
    LUT[ord(_b.lower())] = _i
_REF_OPS = np.array([1, 0, 1, 1, 0, 0, 0, 1, 1, 0], dtype=np.int64)   # M I D N S H P = X B
_Q_OPS = np.array([1, 1, 0, 0, 1, 0, 0, 1, 1, 0], dtype=np.int64)
_M_OPS = np.array([1, 0, 0, 0, 0, 0, 0, 1, 1, 0], dtype=bool)


def fast_strand_counts(bam, contig, L, flt, min_baseq, batch=1000):
    """(4, L) fwd and rev base counts in one pass (vectorised per read from CIGAR blocks);
    same rule as pysam count_coverage: aligned M/=/X bases with quality >= min_baseq."""
    acc = [np.zeros(4 * L, dtype=np.int64), np.zeros(4 * L, dtype=np.int64)]
    keys = [[], []]

    def flush():
        for st in (0, 1):
            if keys[st]:
                acc[st] += np.bincount(np.concatenate(keys[st]), minlength=4 * L)
                keys[st] = []
    n = 0
    for r in bam.fetch(contig):
        if not flt.ok(r) or r.cigartuples is None:
            continue
        seq = r.query_sequence
        if not seq:
            continue
        ops = np.asarray(r.cigartuples, dtype=np.int64)
        op, ln = ops[:, 0], ops[:, 1]
        radv = _REF_OPS[op] * ln
        qadv = _Q_OPS[op] * ln
        rst = r.reference_start + np.cumsum(radv) - radv
        qst = np.cumsum(qadv) - qadv
        m = _M_OPS[op]
        if not m.any():
            continue
        lens, rs, qs = ln[m], rst[m], qst[m]
        tot = int(lens.sum())
        offs = np.arange(tot, dtype=np.int64) - np.repeat(np.cumsum(lens) - lens, lens)
        rpos = np.repeat(rs, lens) + offs
        qpos = np.repeat(qs, lens) + offs
        codes = LUT[np.frombuffer(seq.encode('ascii'), dtype=np.uint8)[qpos]]
        keep = (codes < 4) & (rpos >= 0) & (rpos < L)
        qual = r.query_qualities
        if qual is not None and len(qual):
            keep &= np.frombuffer(qual, dtype=np.uint8)[qpos] >= min_baseq
        keys[1 if r.is_reverse else 0].append(codes[keep] * L + rpos[keep])
        n += 1
        if n % batch == 0:
            flush()
    flush()
    return acc[0].reshape(4, L).astype(np.int32), acc[1].reshape(4, L).astype(np.int32)


def strand_counts(bam, contig, length, flt, min_baseq):
    """(4, L) fwd and rev base counts via pysam count_coverage (slow; reference implementation)."""
    fwd = np.array(bam.count_coverage(contig, 0, length, quality_threshold=min_baseq,
                                      read_callback=lambda r: flt.ok(r) and not r.is_reverse), dtype=np.int32)
    rev = np.array(bam.count_coverage(contig, 0, length, quality_threshold=min_baseq,
                                      read_callback=lambda r: flt.ok(r) and r.is_reverse), dtype=np.int32)
    return fwd, rev


# ---------------------------------------------------------------------- nulls
def encode_contexts(seq, k):
    """Integer code of the k-mer centred on each position (-1 where it contains N or is off the end)."""
    L = len(seq)
    code = np.array([BIDX.get(b, -1) for b in seq], dtype=np.int64)
    h = k // 2
    out = np.zeros(L, dtype=np.int64)
    bad = np.zeros(L, dtype=bool)
    for off in range(-h, h + 1):
        idx = np.arange(L) + off
        valid = (idx >= 0) & (idx < L)
        c = np.full(L, -1, dtype=np.int64)
        c[valid] = code[idx[valid]]
        bad |= c < 0
        out = out * 4 + np.where(c < 0, 0, c)
    out[bad] = -1
    return out


def hp_augment(seq):
    """Per-position code of nearby long homopolymers: run containing the position, the run
    ending just left and the run starting just right; each encoded 0 (< 4) or 1..5 (4..>=8)."""
    L = len(seq)
    run_start = np.zeros(L, dtype=np.int64)
    run_len = np.zeros(L, dtype=np.int64)
    i = 0
    while i < L:
        j = i
        while j + 1 < L and seq[j + 1] == seq[i]:
            j += 1
        run_start[i:j + 1] = i
        run_len[i:j + 1] = j - i + 1
        i = j + 1
    enc = lambda x: np.where(x >= 4, np.minimum(x, 8) - 3, 0)
    own = enc(run_len)
    left = np.zeros(L, dtype=np.int64)
    right = np.zeros(L, dtype=np.int64)
    idx = np.arange(L)
    has_left = run_start > 0
    left[has_left] = enc(run_len[run_start[has_left] - 1])
    nxt = run_start + run_len
    has_right = nxt < L
    right[has_right] = enc(run_len[np.minimum(nxt[has_right], L - 1)])
    return own * 36 + left * 6 + right            # 0..215


HP_CODES = 216


def bb_fit(x, n, trim_p=1e-6, trim_frac=0.02):
    """Moment estimates (mean, rho) of a beta-binomial from counts x of n.

    The most extreme binomial outliers relative to the class mean (P < trim_p) are removed
    before fitting, but at most trim_frac of the class (at least one position). A real
    variant is rare within a context class, so removing a few positions stops it from
    inflating the dispersion and hiding itself; a systematic error shared by many
    positions of the class stays in the fit and widens the null as it should."""
    from scipy.stats import binom
    N = n.sum()
    if N <= 0:
        return 1e-4, 0.0
    p = min(max((x.sum() + 0.5) / (N + 1.0), 1e-6), 0.5)
    pos = np.nonzero(x > 0)[0]
    if len(pos):
        pv = binom.sf(x[pos] - 1, n[pos], p)
        out = pos[pv < trim_p]
        cap = max(1, int(trim_frac * len(x)))
        if len(out) > cap:
            out = pos[np.argsort(pv)[:cap]]
        keep = np.ones(len(x), dtype=bool)
        keep[out] = False
        x, n = x[keep], n[keep]
        N = n.sum()
        if N <= 0:
            return p, 0.0
        p = min(max((x.sum() + 0.5) / (N + 1.0), 1e-6), 0.5)
    s = ((x - n * p) ** 2).sum() / (p * (1 - p))
    den = (n * (n - 1)).sum()
    rho = (s - N) / den if den > 0 else 0.0
    return p, float(min(max(rho, 0.0), 0.5))


def bb_sf(x, n, p, rho):
    """P[X >= x] for X ~ BetaBinomial(n, p, rho) (binomial when rho ~ 0); vectorised."""
    from scipy.stats import betabinom, binom
    x = np.asarray(x)
    n = np.asarray(n)
    p = np.broadcast_to(np.asarray(p, dtype=float), x.shape)
    rho = np.broadcast_to(np.asarray(rho, dtype=float), x.shape)
    out = np.ones(x.shape, dtype=float)
    pos = x > 0
    small = pos & (rho < 1e-7)
    big = pos & ~small
    if small.any():
        out[small] = binom.sf(x[small] - 1, n[small], p[small])
    if big.any():
        a = p[big] * (1 / rho[big] - 1)
        b = (1 - p[big]) * (1 / rho[big] - 1)
        out[big] = betabinom.sf(x[big] - 1, n[big], a, b)
    return np.clip(out, 0.0, 1.0)


class NullModel:
    """Per (strand, context, ALT) beta-binomial parameters with hierarchical fallback.

    Levels: centred 7-, 5-, 3-mer and the base itself, each combined with the homopolymer
    code, then the plain base. A class at a level is used only if it has >= min_positions
    fit positions and >= min_bases bases (the plain-base level always fits). Parameters
    are fitted on non-outlier positions and then applied to every position of the class.
    The error rate is the Jeffreys 97.5 % upper bound of the class rate (conservative for
    classes with few observed errors); the dispersion of classes with < 200 positions is
    floored at the lower quartile of the well-populated classes of the same level."""
    LEVELS = (7, 5, 3, 1, 0)

    def __init__(self, min_positions=30, min_bases=2000):
        self.min_positions, self.min_bases = min_positions, min_bases
        self.summary = []

    def fit(self, data, exclude=None):
        from scipy.stats import beta as beta_dist
        params = [{'p': np.zeros((2, 4, d['L'])), 'rho': np.zeros((2, 4, d['L'])),
                   'level': np.full((2, 4, d['L']), -1, dtype=np.int8)} for d in data]
        for s, key in enumerate(('fwd', 'rev')):
            for alt in range(4):
                assigned = [np.zeros(d['L'], dtype=bool) for d in data]
                for k in self.LEVELS:
                    X, N, C, F, own = [], [], [], [], []
                    for ci, d in enumerate(data):
                        ctx = d[f'ctx{k}']
                        target = d['ok'] & (d['refidx'] != alt) & (ctx >= 0) & ~assigned[ci]
                        idx = np.nonzero(target)[0]
                        if not len(idx):
                            continue
                        X.append(d[key][alt, idx].astype(np.float64))
                        N.append((d['nf'] if s == 0 else d['nr'])[idx].astype(np.float64))
                        C.append(ctx[idx])
                        F.append(np.ones(len(idx), dtype=bool) if exclude is None else ~exclude[ci][s, alt, idx])
                        own.append((ci, idx))
                    if not X:
                        continue
                    X, N, C, F = np.concatenate(X), np.concatenate(N), np.concatenate(C), np.concatenate(F)
                    fitmask = F & (X <= 0.2 * N)
                    order = np.argsort(C, kind='stable')
                    Cs = C[order]
                    bounds = np.flatnonzero(np.diff(Cs)) + 1
                    starts = np.concatenate([[0], bounds])
                    ends = np.concatenate([bounds, [len(Cs)]])
                    res = []
                    for st, en in zip(starts, ends):
                        sel = order[st:en]
                        fm = fitmask[sel]
                        npos, nb = int(fm.sum()), float(N[sel][fm].sum())
                        if k > 0 and (npos < self.min_positions or nb < self.min_bases):
                            res.append(None)
                            continue
                        x_, n_ = X[sel][fm], N[sel][fm]
                        if npos == 0:
                            res.append((1e-3, 0.0, 0, 0.0))
                            continue
                        p_, rho_ = bb_fit(x_, n_)
                        # conservative rate: Jeffreys 97.5 % upper bound on the pooled class rate
                        xs_, ns_ = float(x_.sum()), float(n_.sum())
                        p_up = float(beta_dist.ppf(0.975, xs_ + 0.5, ns_ - xs_ + 0.5))
                        res.append((max(p_up, 1e-6), rho_, npos, xs_ / max(1.0, ns_)))
                    big = [r[1] for r in res if r is not None and r[2] >= 200]
                    rho_floor = float(np.quantile(big, 0.25)) if big else 0.0
                    Pv = np.full(len(C), np.nan)
                    Rv = np.full(len(C), np.nan)
                    for (st, en), r in zip(zip(starts, ends), res):
                        if r is None:
                            continue
                        p_, rho_, npos, raw = r
                        if npos < 200:
                            rho_ = max(rho_, rho_floor)
                        sel = order[st:en]
                        Pv[sel] = p_
                        Rv[sel] = rho_
                        self.summary.append((s, alt, k, int(Cs[st]), npos, raw, p_, rho_))
                    off = 0
                    for ci, idx in own:
                        m = len(idx)
                        pp, rr = Pv[off:off + m], Rv[off:off + m]
                        got = ~np.isnan(pp)
                        params[ci]['p'][s, alt, idx[got]] = pp[got]
                        params[ci]['rho'][s, alt, idx[got]] = rr[got]
                        params[ci]['level'][s, alt, idx[got]] = k
                        assigned[ci][idx[got]] = True
                        off += m
        return params


# ---------------------------------------------------------------------- calls
def bh(p):
    p = np.asarray(p, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


def hp_len(seq, i):
    b = seq[i]
    s = e = i
    while s > 0 and seq[s - 1] == b:
        s -= 1
    while e + 1 < len(seq) and seq[e + 1] == b:
        e += 1
    return e - s + 1


def site_alleles(bam, contig, pos0, flt, min_baseq, max_depth=1_000_000):
    """read name -> base at pos0 for filtered reads (pileup)."""
    out = {}
    for col in bam.pileup(contig, pos0, pos0 + 1, truncate=True, stepper='nofilter', ignore_orphans=False,
                          min_base_quality=0, max_depth=max_depth, ignore_overlaps=False):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            r = pr.alignment
            if not flt.ok(r):
                continue
            q = r.query_qualities
            if q is not None and q[pr.query_position] < min_baseq:
                continue
            out[r.query_name] = r.query_sequence[pr.query_position]
    return out


def fisher_greater(a, b, c, d):
    from scipy.stats import hypergeom
    # P[X >= a] with X ~ Hypergeom(M=a+b+c+d, n=a+b (alt at i), N=a+c (alt at j))
    return float(hypergeom.sf(a - 1, a + b + c + d, a + b, a + c))


def link_groups(cands, alleles, a):
    """Union-find over linked candidate pairs; returns list of groups (lists of candidate indices)."""
    parent = list(range(len(cands)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    by_contig = {}
    for i, c in enumerate(cands):
        by_contig.setdefault(c['contig'], []).append(i)
    tested = linked = 0
    A = [{n for n, b in alleles[i].items() if b == cands[i]['alt']} for i in range(len(cands))]
    R = [{n for n, b in alleles[i].items() if b == cands[i]['ref']} for i in range(len(cands))]
    for ctg, idx in by_contig.items():
        idx.sort(key=lambda i: cands[i]['pos'])
        for x in range(len(idx)):
            i = idx[x]
            Ai, Ri = A[i], R[i]
            if len(Ai) < 3:
                continue
            for y in range(x + 1, len(idx)):
                j = idx[y]
                if cands[j]['pos'] - cands[i]['pos'] > a.link_max_dist:
                    break
                Aj, Rj = A[j], R[j]
                aa, ab, ba, bb = len(Ai & Aj), len(Ai & Rj), len(Ri & Aj), len(Ri & Rj)
                if aa + ab + ba + bb < 10 or aa < 3:
                    continue
                tested += 1
                orr = (aa + 0.5) * (bb + 0.5) / ((ab + 0.5) * (ba + 0.5))
                if orr < a.link_min_or:
                    continue
                if fisher_greater(aa, ab, ba, bb) <= a.link_p:
                    linked += 1
                    parent[find(i)] = find(j)
    groups = {}
    for i in range(len(cands)):
        groups.setdefault(find(i), []).append(i)
    log(f'linkage: {tested} pairs tested, {linked} linked')
    return [g for g in groups.values() if len(g) >= 2]


def group_reads(group, cands, alleles):
    votes = {}
    for i in group:
        for n, b in alleles[i].items():
            v = votes.setdefault(n, [0, 0])
            if b == cands[i]['alt']:
                v[0] += 1
            elif b == cands[i]['ref']:
                v[1] += 1
    minor = major = amb = 0
    for alt_n, ref_n in votes.values():
        cov = alt_n + ref_n
        if cov < 2:
            continue
        if alt_n >= 2 and alt_n * 3 >= 2 * cov:
            minor += 1
        elif ref_n >= 2 and ref_n * 3 >= 2 * cov:
            major += 1
        else:
            amb += 1
    return minor, major, amb


# ----------------------------------------------------------------------- main
def cmd_call(a):
    import pysam
    t0 = time.time()
    os.makedirs(a.outdir, exist_ok=True)
    recs = read_fasta(a.ref)
    bam = pysam.AlignmentFile(a.bam)
    names = set(bam.references)
    recs = [r for r in recs if r[0] in names]
    log(f'{a.sample}: {len(recs)} contigs, {sum(len(s) for _, s, _ in recs)} bp')
    flt = ReadFilter(a.min_mapq, 1.0, a.seed)
    md = None
    if a.max_depth > 0:
        md = mean_depth(bam, recs, flt)
        flt.keep_frac = min(1.0, a.max_depth / md) if md > 0 else 1.0
        log(f'mean depth {md:.0f}x; keeping fraction {flt.keep_frac:.4f} of reads')
    if a.repeats:
        reps = load_bed(a.repeats)
    else:
        reps = repeat_intervals(recs)
    feats = load_gff(a.gff, gff_name_map(recs, a.gff_fasta))
    if a.gff and not any(n in feats for n, _, _ in recs):
        log('WARNING: no GFF feature is on a consensus contig; effects will be missing (pass --gff-fasta?)')
    data = []
    for name, seq, hdr in recs:
        L = len(seq)
        fwd, rev = fast_strand_counts(bam, name, L, flt, a.min_baseq)
        refidx = np.array([BIDX.get(b, -1) for b in seq], dtype=np.int8)
        mask = np.zeros(L, dtype=bool)
        for s, e in reps.get(name, []):
            mask[s:e] = True
        mask[:a.end_mask] = True
        mask[max(0, L - a.end_mask):] = True
        mask |= refidx < 0
        nf, nr = fwd.sum(axis=0), rev.sum(axis=0)
        ok = ~mask & (nf >= a.min_strand_depth) & (nr >= a.min_strand_depth)
        d = dict(name=name, seq=seq, L=L, fwd=fwd, rev=rev, refidx=refidx, mask=mask, ok=ok, nf=nf, nr=nr)
        aug = hp_augment(seq)
        for k in NullModel.LEVELS:
            if k == 0:
                d['ctx0'] = refidx.astype(np.int64)
                continue
            base = encode_contexts(seq, k) if k > 1 else refidx.astype(np.int64)
            d[f'ctx{k}'] = np.where(base >= 0, base * HP_CODES + aug, -1)
        data.append(d)
        log(f'{name}: counted; assessed {int(ok.sum())} of {L} positions (masked {int(mask.sum())})')
    # ---- null model (capped outlier trimming inside each class fit)
    nm = NullModel(a.min_class_positions, a.min_class_bases)
    params = nm.fit(data)
    log('null model fitted')
    # ---- tests (vectorised): cells with >= min ALT reads on both strands
    T_ci, T_pos, T_alt, T_p = [], [], [], []
    m_tests = 0
    for ci, (d, pr) in enumerate(zip(data, params)):
        m_tests += 3 * int(d['ok'].sum())
        for alt in range(4):
            xf, xr = d['fwd'][alt], d['rev'][alt]
            cand = d['ok'] & (d['refidx'] != alt) & (xf >= a.min_alt_per_strand) & (xr >= a.min_alt_per_strand)
            idx = np.nonzero(cand)[0]
            if not len(idx):
                continue
            pf = bb_sf(xf[idx], d['nf'][idx], pr['p'][0, alt, idx], pr['rho'][0, alt, idx])
            pr_ = bb_sf(xr[idx], d['nr'][idx], pr['p'][1, alt, idx], pr['rho'][1, alt, idx])
            T_ci.append(np.full(len(idx), ci)); T_pos.append(idx); T_alt.append(np.full(len(idx), alt))
            T_p.append(np.maximum(pf, pr_))
    if T_p:
        T_ci, T_pos, T_alt, pvals = (np.concatenate(x) for x in (T_ci, T_pos, T_alt, T_p))
    else:
        T_ci = T_pos = T_alt = np.zeros(0, dtype=np.int64); pvals = np.zeros(0)
    # Benjamini-Hochberg over the full family; untested cells have p = 1
    q = np.ones(len(pvals))
    if len(pvals):
        allp = np.concatenate([pvals, np.ones(max(0, m_tests - len(pvals)))])
        q = bh(allp)[:len(pvals)]
    called = q <= a.fdr
    rows = list(zip(T_ci.tolist(), T_pos.tolist(), T_alt.tolist()))
    # strand bias: reject a call whose per-strand fractions differ > --sb-ratio-fold AND Fisher p < --sb-p
    n_sb = 0
    if called.any():
        from scipy.stats import fisher_exact
        for k in np.nonzero(called)[0].tolist():
            ci, i, alt = rows[k]
            d = data[ci]
            xf, xr, nf, nr = int(d['fwd'][alt, i]), int(d['rev'][alt, i]), int(d['nf'][i]), int(d['nr'][i])
            ff, fr = xf / max(1, nf), xr / max(1, nr)
            ratio = min(ff, fr) / max(ff, fr) if max(ff, fr) > 0 else 1.0
            if ratio < 1.0 / a.sb_ratio_fold:
                if fisher_exact([[xf, nf - xf], [xr, nr - xr]])[1] < a.sb_p:
                    called[k] = False
                    n_sb += 1
    log(f'strand-bias filter removed {n_sb} calls')
    # context enrichment: systematic errors recur across positions sharing a context; real
    # variants are scattered. Calls whose (5-mer + homopolymer code, ALT) class holds >= 3 calls
    # and more than chance allows (Poisson, Bonferroni over classes) are flagged 'systematic'.
    systematic = np.zeros(len(called), dtype=bool)
    if called.any():
        from scipy.stats import poisson
        key_pos = {}
        for ci, d in enumerate(data):
            kk = d['ctx5'][d['ok']]
            for alt in range(4):
                sel = d['refidx'][d['ok']] != alt
                u, cnt = np.unique(kk[sel] * 4 + alt, return_counts=True)
                for uu, cc in zip(u.tolist(), cnt.tolist()):
                    key_pos[uu] = key_pos.get(uu, 0) + cc
        tot_cells = sum(key_pos.values())
        ck = {}
        for k in np.nonzero(called)[0].tolist():
            ci, i, alt = rows[k]
            ck.setdefault(int(data[ci]['ctx5'][i]) * 4 + alt, []).append(k)
        n_called = int(called.sum())
        thr = a.ctx_p / max(1, len(key_pos))
        for key, ks in ck.items():
            if len(ks) < 3:
                continue
            lam = n_called * key_pos.get(key, 1) / max(1, tot_cells)
            if poisson.sf(len(ks) - 1, lam) <= thr:
                systematic[ks] = True
        log(f'context-enrichment filter flagged {int(systematic.sum())} calls in '
            f'{sum(1 for ks in ck.values() if systematic[ks].any())} context classes')
    called &= ~systematic
    log(f'{m_tests} tests; {len(pvals)} cells with >= {a.min_alt_per_strand} ALT reads per strand; '
        f'{int(called.sum())} called at q <= {a.fdr}')
    # ---- candidates for linkage
    cands = []
    for k in np.nonzero(called | (pvals <= a.candidate_p))[0].tolist():
        ci, i, alt = rows[k]
        d = data[ci]
        cands.append(dict(k=k, ci=ci, contig=d['name'], pos=i, ref=d['seq'][i], alt=BASES[alt],
                          called=bool(called[k])))
    # also sub-threshold candidates with >= 2 ALT reads per strand
    for ci, (d, pr) in enumerate(zip(data, params)):
        for alt in range(4):
            xf, xr = d['fwd'][alt], d['rev'][alt]
            cand = d['ok'] & (d['refidx'] != alt) & (xf >= 2) & (xr >= 2) & \
                   ((xf < a.min_alt_per_strand) | (xr < a.min_alt_per_strand))
            idx = np.nonzero(cand)[0]
            if not len(idx):
                continue
            pf = bb_sf(xf[idx], d['nf'][idx], pr['p'][0, alt, idx], pr['rho'][0, alt, idx])
            prr = bb_sf(xr[idx], d['nr'][idx], pr['p'][1, alt, idx], pr['rho'][1, alt, idx])
            for j in np.nonzero(np.maximum(pf, prr) <= a.candidate_p)[0]:
                cands.append(dict(k=None, ci=ci, contig=d['name'], pos=int(idx[j]), ref=d['seq'][idx[j]],
                                  alt=BASES[alt], called=False, p=float(max(pf[j], prr[j]))))
    if len(cands) > a.max_link_candidates:
        log(f'{len(cands)} linkage candidates > {a.max_link_candidates}; keeping called sites and the lowest p-values')
        cands.sort(key=lambda c: (not c['called'], c.get('p', pvals[c['k']] if c['k'] is not None else 1)))
        cands = cands[:a.max_link_candidates]
    alleles = [site_alleles(bam, c['contig'], c['pos'], flt, a.min_baseq) for c in cands]
    groups = link_groups(cands, alleles, a) if len(cands) >= 2 else []
    site_group, group_rows = {}, []
    gid = 0
    for g in sorted(groups, key=lambda g: (cands[g[0]]['contig'], min(cands[i]['pos'] for i in g))):
        if not any(cands[i]['called'] for i in g):
            continue
        gid += 1
        name = f'G{gid}'
        minor, major, amb = group_reads(g, cands, alleles)
        lo, hi = wilson(minor, minor + major)
        poss = sorted(cands[i]['pos'] for i in g)
        group_rows.append([name, cands[g[0]]['contig'], poss[0] + 1, poss[-1] + 1, len(g),
                           sum(cands[i]['called'] for i in g), minor, major, amb,
                           f'{minor / max(1, minor + major):.5f}', f'{lo:.5f}', f'{hi:.5f}',
                           ','.join(str(p + 1) for p in poss)])
        for i in g:
            site_group[(cands[i]['contig'], cands[i]['pos'], cands[i]['alt'])] = name
    # ---- site table
    out_sites = []
    keep_k = set(np.nonzero(called)[0].tolist())
    for c in cands:
        if c['k'] is not None and (c['contig'], c['pos'], c['alt']) in site_group:
            keep_k.add(c['k'])
    for k in sorted(keep_k):
        ci, i, alt = rows[k]
        out_sites.append((k, ci, i, alt, 'called' if called[k] else 'linked'))
    for c in cands:
        key = (c['contig'], c['pos'], c['alt'])
        if c['k'] is None and key in site_group:
            out_sites.append((None, c['ci'], c['pos'], BASES.index(c['alt']), 'linked'))
    cols = ['contig', 'pos', 'ref', 'alt', 'status', 'fraction', 'ci_low', 'ci_high', 'raw_fraction', 'background',
            'alt_fwd', 'alt_rev', 'depth_fwd', 'depth_rev', 'p', 'q', 'context', 'hp_len', 'group', 'locus_tag',
            'gene', 'product', 'effect', 'aa_change']
    lines = []
    for k, ci, i, alt, status in sorted(out_sites, key=lambda x: (x[1], x[2], x[3])):
        d = data[ci]
        pr = params[ci]
        xf, xr, nf, nr = int(d['fwd'][alt, i]), int(d['rev'][alt, i]), int(d['nf'][i]), int(d['nr'][i])
        n, x = nf + nr, xf + xr
        e = (pr['p'][0, alt, i] * nf + pr['p'][1, alt, i] * nr) / max(1, n)
        raw = x / max(1, n)
        lo, hi = wilson(x, n)
        adj = lambda v: min(1.0, max(0.0, (v - e) / (1 - e)))
        if raw > 0.5:
            status = 'majority_alt'
        pv = pvals[k] if k is not None else float('nan')
        qv = q[k] if k is not None else float('nan')
        loc, gene, prod, eff, hg = snv_effect(d['seq'], feats.get(d['name'], []), i + 1, d['seq'][i], BASES[alt])
        ctx = d['seq'][max(0, i - 3):i + 4]
        lines.append([d['name'], i + 1, d['seq'][i], BASES[alt], status, f'{adj(raw):.5f}', f'{adj(lo):.5f}',
                      f'{adj(hi):.5f}', f'{raw:.5f}', f'{e:.6f}', xf, xr, nf, nr, f'{pv:.3g}', f'{qv:.3g}', ctx,
                      hp_len(d['seq'], i), site_group.get((d['name'], i, BASES[alt]), '-'), loc, gene, prod, eff, hg])
    with open(os.path.join(a.outdir, 'sites.tsv'), 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in lines:
            fh.write('\t'.join(str(x).replace('\t', ' ') for x in r) + '\n')
    with open(os.path.join(a.outdir, 'groups.tsv'), 'w') as fh:
        fh.write('group\tcontig\tstart\tend\tn_sites\tn_called\treads_minor\treads_major\treads_ambiguous\t'
                 'fraction\tci_low\tci_high\tpositions\tgenes\n')
        gg = group_genes(lines)
        for r in group_rows:
            fh.write('\t'.join(map(str, r)) + '\t' + gg.get(str(r[0]), '') + '\n')
    # ---- known versions
    known = known_versions(a, data, bam, flt) if a.version else []
    # ---- detection floor: smallest fraction that reaches p <= fdr / m at a site's depth and null
    floors = []
    rng = np.random.default_rng(1)
    for d, pr in zip(data, params):
        idx = np.nonzero(d['ok'])[0]
        if not len(idx):
            continue
        pick = rng.choice(idx, size=min(2000, len(idx)), replace=False)
        thr = a.fdr / max(1, m_tests)
        for i in pick:
            alt = (int(d['refidx'][i]) + 1) % 4
            nf, nr = int(d['nf'][i]), int(d['nr'][i])
            for f in (0.001, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5):
                xf, xr = max(a.min_alt_per_strand, round(f * nf)), max(a.min_alt_per_strand, round(f * nr))
                if xf > nf or xr > nr:
                    continue
                pv = max(bb_sf(np.array([xf]), np.array([nf]), pr['p'][0, alt, i], pr['rho'][0, alt, i])[0],
                         bb_sf(np.array([xr]), np.array([nr]), pr['p'][1, alt, i], pr['rho'][1, alt, i])[0])
                if pv <= thr:
                    floors.append(f)
                    break
            else:
                floors.append(1.0)
    floors = np.array(floors) if floors else np.array([1.0])
    depth_all = np.concatenate([(d['nf'] + d['nr'])[d['ok']] for d in data]) if data else np.array([0])
    summ = dict(sample=a.sample, positions_total=sum(d['L'] for d in data),
                positions_assessed=int(sum(d['ok'].sum() for d in data)),
                positions_masked=int(sum(d['mask'].sum() for d in data)),
                positions_low_depth=int(sum((~d['mask'] & ~d['ok']).sum() for d in data)),
                median_depth=int(np.median(depth_all)) if len(depth_all) else 0,
                n_sites=len(lines), n_called=sum(1 for r in lines if r[4] == 'called'),
                n_linked=sum(1 for r in lines if r[4] == 'linked'),
                n_majority=sum(1 for r in lines if r[4] == 'majority_alt'), n_groups=len(group_rows),
                floor_p10=float(np.quantile(floors, 0.1)), floor_median=float(np.median(floors)),
                floor_p90=float(np.quantile(floors, 0.9)), tests=m_tests, fdr=a.fdr,
                read_fraction_used=flt.keep_frac, mean_depth_before_cap=md if md is not None else '',
                known_versions=len(known), engine_version=VERSION, seconds=int(time.time() - t0))
    with open(os.path.join(a.outdir, 'summary.tsv'), 'w') as fh:
        fh.write('\t'.join(summ) + '\n' + '\t'.join(str(v) for v in summ.values()) + '\n')
    # null model summary
    lv = {}
    for s_, alt, k, code, npos, raw, p_, rho_ in nm.summary:
        lv.setdefault(k, []).append((raw, p_, rho_, npos))
    with open(os.path.join(a.outdir, 'null_model.tsv'), 'w') as fh:
        fh.write('level\tclasses\tpositions\traw_rate_median\tused_rate_p10\tused_rate_median\tused_rate_p90\trho_median\trho_p90\n')
        for k in sorted(lv, reverse=True):
            v = np.array(lv[k])
            fh.write(f'{k}\t{len(v)}\t{int(v[:, 3].sum())}\t{np.median(v[:, 0]):.3g}\t{np.quantile(v[:, 1], .1):.3g}\t'
                     f'{np.median(v[:, 1]):.3g}\t{np.quantile(v[:, 1], .9):.3g}\t{np.median(v[:, 2]):.3g}\t'
                     f'{np.quantile(v[:, 2], .9):.3g}\n')
    log(f'done: {summ["n_sites"]} sites ({summ["n_called"]} called, {summ["n_linked"]} linked, '
        f'{summ["n_majority"]} majority-ALT), {summ["n_groups"]} groups; floor median {summ["floor_median"]}')


def known_versions(a, data, bam, flt):
    """Genotype SNP differences between each supplied genome version and the consensus."""
    import mappy
    rows, vrows = [], []
    for vpath in a.version:
        vrecs = read_fasta(vpath)
        vname = os.path.splitext(os.path.basename(vpath))[0]
        sites = []
        for d in data:
            al = mappy.Aligner(seq=d['seq'], preset='asm5')
            for vn, vs, _ in vrecs:
                for h in al.map(vs, cs=True):
                    if not h.is_primary:
                        continue
                    r = h.r_st
                    for m in re.finditer(r'(:\d+|\*[a-z][a-z]|[+-][a-z]+)', h.cs):
                        t = m.group(1)
                        if t[0] == ':':
                            r += int(t[1:])
                        elif t[0] == '*':
                            alt = t[2].upper() if h.strand == 1 else t[2].upper()
                            if not d['mask'][r]:
                                sites.append((d, r, alt))
                            r += 1
                        elif t[0] == '-':
                            r += len(t) - 1
        if not sites:
            vrows.append([vname, 0, 0, 0, 0, '', '', '', 'no SNP differences from the consensus outside masked sequence'])
            continue
        votes = {}
        for d, r, alt in sites:
            al = site_alleles(bam, d['name'], r, flt, a.min_baseq)
            nv = sum(1 for b in al.values() if b == alt)
            nc = sum(1 for b in al.values() if b == d['seq'][r])
            lo, hi = wilson(nv, nv + nc)
            rows.append([vname, d['name'], r + 1, d['seq'][r], alt, nv, nc, f'{nv / max(1, nv + nc):.5f}',
                         f'{lo:.5f}', f'{hi:.5f}'])
            for n, b in al.items():
                v = votes.setdefault(n, [0, 0])
                v[0] += b == alt
                v[1] += b == d['seq'][r]
        rv = sum(1 for x, y in votes.values() if x > y and x >= 1)
        rc_ = sum(1 for x, y in votes.values() if y > x and y >= 1)
        multi = sum(1 for x, y in votes.values() if x + y >= 2)
        lo, hi = wilson(rv, rv + rc_)
        vrows.append([vname, len(sites), rv, rc_, multi, f'{rv / max(1, rv + rc_):.5f}', f'{lo:.5f}', f'{hi:.5f}', ''])
    with open(os.path.join(a.outdir, 'known_sites.tsv'), 'w') as fh:
        fh.write('version\tcontig\tpos\tconsensus\tversion_allele\treads_version\treads_consensus\tfraction\tci_low\tci_high\n')
        for r in rows:
            fh.write('\t'.join(map(str, r)) + '\n')
    with open(os.path.join(a.outdir, 'known_versions.tsv'), 'w') as fh:
        fh.write('version\tsnp_sites\treads_version\treads_consensus\treads_multi_site\tfraction\tci_low\tci_high\tnote\n')
        for r in vrows:
            fh.write('\t'.join(map(str, r)) + '\n')
    return vrows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version-info', action='version', version=VERSION)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('call')
    p.add_argument('--bam', required=True); p.add_argument('--ref', required=True)
    p.add_argument('--sample', required=True); p.add_argument('--outdir', required=True)
    p.add_argument('--gff', default=''); p.add_argument('--gff-fasta', default='', help='annotation FASTA (maps renamed contigs)')
    p.add_argument('--repeats', default='')
    p.add_argument('--version', action='append', default=[], help='genome version FASTA (repeatable)')
    p.add_argument('--min-baseq', type=int, default=10); p.add_argument('--min-mapq', type=int, default=20)
    p.add_argument('--min-strand-depth', type=int, default=20)
    p.add_argument('--min-alt-per-strand', type=int, default=3)
    p.add_argument('--fdr', type=float, default=0.01); p.add_argument('--max-depth', type=float, default=0)
    p.add_argument('--end-mask', type=int, default=500); p.add_argument('--seed', default='minor1')
    p.add_argument('--candidate-p', type=float, default=1e-4)
    p.add_argument('--link-max-dist', type=int, default=50000); p.add_argument('--link-min-or', type=float, default=20)
    p.add_argument('--link-p', type=float, default=1e-8)
    p.add_argument('--max-link-candidates', type=int, default=20000)
    p.add_argument('--min-class-positions', type=int, default=30)
    p.add_argument('--min-class-bases', type=int, default=2000)
    p.add_argument('--sb-ratio-fold', type=float, default=3.0)
    p.add_argument('--sb-p', type=float, default=1e-3)
    p.add_argument('--ctx-p', type=float, default=1e-3, help='family-wise p for context enrichment')
    a = ap.parse_args()
    if a.cmd == 'call':
        cmd_call(a)


if __name__ == '__main__':
    main()
