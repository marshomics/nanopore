#!/usr/bin/env python3
"""
sv_genotype.py - structural heterogeneity in prokaryotic long-read data.

Every candidate locus is represented by explicit allele sequences: a shared
left flank, one version of the variable segment per allele, and a shared
right flank. Candidates come from (i) differences between supplied genome
versions of the organism, (ii) Sniffles2 calls on each sample's reads, and
(iii) size-matched fake events at random unique positions (empirical null).

Reads at a locus are re-aligned to every allele with mappy. A read counts for
an allele only if a single alignment covers that allele's window without an
internal indel as large as the allele difference. Windows have the same
length for every allele (the shorter allele's window is padded with flank),
so a genome copy carrying either allele has the same chance of producing a
counted read. Anchors at both ends of the window contain >= --anchor bp of
sequence that is unique in the genome and within every allele, so repeats
and tandem duplications cannot place a read on the wrong allele.

Fractions are reported with Wilson 95% intervals. A per-sample detection
floor is the highest alt fraction seen at the null loci (at least
--min-fraction); an allele is called present when it has >= --min-alt-reads
reads and its lower 95% bound exceeds the floor.

Events whose allele windows are longer than the reads (large deletions,
tandem amplifications, excisable elements, long inversions) are genotyped at
their breakpoint junctions instead (--junction-mode auto; 'always' runs both
methods). For a segment, four junction alleles are built from the reported,
read-supported breakpoints: left flank|segment start, segment end|right
flank, segment end|segment start (tandem copy or circular form) and left
flank|right flank (segment absent). All four share one window geometry that
reaches a unique anchor on both sides of any repeat next to the junction.
Each region of a read that spans a junction window supports one junction
allele. Estimates: fraction of genomes lacking the segment, and end-to-start
junctions per genome (tandem copies or circles, not distinguishable by
junctions alone), each with a Wilson interval and a junction-specific null
floor. Inversions use reference/inverted junctions at both ends.
"""
import argparse
import hashlib
import json
import math
import os
import random
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import mappy
import pysam

VERSION = "1.4"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
NULL_BINS = (500, 1000, 2000, 5000, 10000, 20000, 50000)


def rc(s):
    return s.translate(COMP)[::-1]


def log(msg):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), msg, flush=True)


def wilson(k, n, z=1.96):
    if n <= 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, (c - h) / d), min(1.0, (c + h) / d)


def name_hash(name):
    return int(hashlib.md5(name.encode()).hexdigest()[:12], 16)


# ----------------------------------------------------------------- sequences
class Contig:
    def __init__(self, name, comment, seq):
        self.name = name
        self.comment = comment or ""
        self.seq = seq.upper()
        self.len = len(self.seq)
        c = self.comment.lower()
        self.circular = "circular=true" in c or "topology=circular" in c
        self._tri = None

    def sub(self, start, end):
        """0-based half-open slice; wraps around the ends of circular contigs.
        Returns None for linear contigs when the window leaves the contig."""
        if end - start > self.len:
            return None
        if self.circular:
            if self._tri is None:
                self._tri = self.seq * 3
            if start < -self.len or end > 2 * self.len:
                return None
            return self._tri[start + self.len:end + self.len]
        if start < 0 or end > self.len:
            return None
        return self.seq[start:end]


def read_fasta(path):
    out = []
    with pysam.FastxFile(path) as fh:
        for e in fh:
            out.append(Contig(e.name, e.comment, e.sequence))
    return out


def merge_intervals(ivs):
    out = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def self_repeats(seq, chunk=1000, step=500, min_len=300, min_ident=0.85):
    """Intervals of `seq` that align elsewhere in `seq` (either strand)."""
    if len(seq) < chunk:
        return []
    al = mappy.Aligner(seq=seq, preset="asm20", best_n=50)
    if not al:
        return []
    ivs = []
    last = len(seq) - chunk
    starts = list(range(0, last + 1, step))
    if starts[-1] != last:
        starts.append(last)
    for st in starts:
        q = seq[st:st + chunk]
        for h in al.map(q):
            if h.strand == 1 and abs(h.r_st - (st + h.q_st)) <= 50:
                continue
            if h.blen < min_len or h.mlen < min_ident * h.blen:
                continue
            ivs.append((st + h.q_st, st + h.q_en))
    return merge_intervals(ivs)


class RepeatMap:
    def __init__(self, contigs):
        self.ivs = {}
        for c in contigs:
            if c.circular and c.len > 4000:
                raw = self_repeats(c.seq + c.seq[:2000])   # include the origin junction
                self.ivs[c.name] = merge_intervals([(s, min(e, c.len)) for s, e in raw if s < c.len])
            else:
                self.ivs[c.name] = self_repeats(c.seq)

    def masked_fraction(self, contig, start, end):
        tot = 0
        for s, e in self.ivs.get(contig, []):
            if e <= start:
                continue
            if s >= end:
                break
            tot += min(e, end) - max(s, start)
        return tot / max(1, end - start)

    def mask_array(self, contig, start, end, clen, circular):
        """Boolean list for positions start..end-1 (may wrap on circular contigs)."""
        out = [False] * (end - start)
        for s, e in self.ivs.get(contig, []):
            for shift in ((-clen, 0, clen) if circular else (0,)):
                lo, hi = max(s + shift, start), min(e + shift, end)
                for p in range(lo, hi):
                    out[p - start] = True
        return out


# ----------------------------------------------------------------- loci
class Allele:
    def __init__(self, name, seq, vs, ve, label=""):
        self.name, self.seq, self.vs, self.ve = name, seq, vs, ve
        self.label = label
        self.win = None

    @property
    def vlen(self):
        return self.ve - self.vs


class Locus:
    def __init__(self, lid, source, kind, alleles, desc, ext_mask_left=None, ext_mask_right=None):
        self.id, self.source, self.kind, self.alleles, self.desc = lid, source, kind, alleles, desc
        self.FL = alleles[0].vs
        self.FR = len(alleles[0].seq) - alleles[0].ve
        for a in alleles:
            assert a.vs == self.FL and len(a.seq) - a.ve == self.FR, "alleles must share flank lengths"
        self.ext_mask_left = ext_mask_left
        self.ext_mask_right = ext_mask_right
        self.A_left = self.A_right = None
        self.W = None
        self.status_note = ""
        self.repeat_limited = False
        self.max_indel = None
        self.meta = {}

    @property
    def lmax(self):
        return max(a.vlen for a in self.alleles)

    @property
    def dlen(self):
        ls = [a.vlen for a in self.alleles]
        return max(ls) - min(ls)


def compute_windows(locus, anchor, margin=300):
    """Repeat-aware anchors and equal-length windows for every allele."""
    FL, FR = locus.FL, locus.FR
    lm = list(locus.ext_mask_left) if locus.ext_mask_left else [False] * FL
    rmk = list(locus.ext_mask_right) if locus.ext_mask_right else [False] * FR
    for a in locus.alleles:
        for s, e in self_repeats(a.seq):
            for p in range(max(s, 0), min(e, FL)):
                lm[p] = True
            for p in range(max(s, a.ve), min(e, len(a.seq))):
                rmk[p - a.ve] = True
    A_left = A_right = None
    run = 0
    for d in range(1, FL + 1):
        run = 0 if lm[FL - d] else run + 1
        if run >= anchor:
            A_left = d
            break
    run = 0
    for d in range(1, FR + 1):
        run = 0 if rmk[d - 1] else run + 1
        if run >= anchor:
            A_right = d
            break
    if A_left is None or A_right is None:
        locus.repeat_limited = True
        locus.status_note = "no unique anchor within flanks"
        return False
    locus.A_left, locus.A_right = A_left, A_right
    lmax = locus.lmax
    for a in locus.alleles:
        pad = lmax - a.vlen
        pl = pad // 2
        pr = pad - pl
        ws, we = a.vs - A_left - pl, a.ve + A_right + pr
        if ws < margin or we > len(a.seq) - margin:
            locus.repeat_limited = True
            locus.status_note = "flanks too short for anchors plus padding"
            return False
        a.win = (ws, we)
    locus.W = A_left + A_right + lmax
    d = locus.dlen
    locus.max_indel = 500 if d == 0 else max(50, min(500, d // 2))
    return True


def make_locus_from_contig(contig, rmap, lid, source, kind, start, end, alt_segment, flank, ref_label, alt_label, desc):
    """Two-allele locus from a contig: REF = contig[start:end], ALT = alt_segment."""
    left = contig.sub(start - flank, start)
    right = contig.sub(end, end + flank)
    mid = contig.sub(start, end)
    if left is None or right is None or mid is None:
        return None
    ref = Allele("REF", left + mid + right, flank, flank + len(mid), ref_label)
    alt = Allele("ALT", left + alt_segment + right, flank, flank + len(alt_segment), alt_label)
    ml = rmap.mask_array(contig.name, start - flank, start, contig.len, contig.circular)
    mr = rmap.mask_array(contig.name, end, end + flank, contig.len, contig.circular)
    L = Locus(lid, source, kind, [ref, alt], desc, ml, mr)
    L.meta.update(contig=contig.name, start=start, end=end)
    return L


def locus_from_full(contig, rmap, lid, source, kind, x0, x1, named, desc):
    """Locus from full allele sequences over contig[x0:x1]; shared flanks = common prefix / suffix."""
    seqs = [sq for _, sq, _ in named]
    if any(sq is None for sq in seqs):
        return None
    FL = len(os.path.commonprefix(seqs))
    FR = len(os.path.commonprefix([sq[::-1] for sq in seqs]))
    FR = min(FR, min(len(sq) for sq in seqs) - FL)
    alleles = [Allele(nm, sq, FL, len(sq) - FR, lab) for nm, sq, lab in named]
    ml = rmap.mask_array(contig.name, x0, x0 + FL, contig.len, contig.circular)
    mr = rmap.mask_array(contig.name, x1 - FR, x1, contig.len, contig.circular)
    L = Locus(lid, source, kind, alleles, desc, ml, mr)
    L.meta.update(contig=contig.name, start=x0 + FL, end=x1 - FR)
    return L


def direct_repeat_pair(contig, s, e, X=6000, min_len=100, min_ident=0.9, tol=150):
    """Direct repeat copies D1 at breakpoint s and D2 at breakpoint e: (a1, a2, b1, b2) or None.

    Split alignments in a repeat may overlap in read coordinates (both parts claim the repeat), so a
    caller can report breakpoints at the outer (or inner) repeat edges although the molecules carry
    one copy at the junction. The caller then adds the one-copy allele as an alternative."""
    X = min(X, (e - s) // 2)
    if X < min_len:
        return None
    A, B = contig.sub(s - X, s + X), contig.sub(e - X, e + X)
    if A is None or B is None:
        return None
    best = None
    for h in mappy.Aligner(seq=B, preset="asm10").map(A):
        if h.strand != 1 or h.q_en - h.q_st < min_len or h.mlen < min_ident * h.blen:
            continue
        a1, a2, b1, b2 = s - X + h.q_st, s - X + h.q_en, e - X + h.r_st, e - X + h.r_en
        if a1 - tol <= s <= a2 + tol and b1 - tol <= e <= b2 + tol and a2 <= b1:
            if best is None or a2 - a1 > best[1] - best[0]:
                best = (a1, a2, b1, b2)
    return best


# ------------------------------------------------- known genome versions
def normalise_to(base, other, probe=3000):
    """Rotate / reverse-complement `other` so that it starts where `base` starts."""
    for o, strand in ((other, "+"), (rc(other), "-")):
        al = mappy.Aligner(seq=o, preset="asm5", best_n=5)
        for off in (0, len(base) // 3, 2 * len(base) // 3):
            q = base[off:off + probe]
            hits = [h for h in al.map(q) if h.q_en - h.q_st >= 0.9 * len(q)]
            if len(hits) == 1 and hits[0].strand == 1:
                p = (hits[0].r_st - hits[0].q_st - off) % len(o)
                return o[p:] + o[:p], strand, p
    return None, None, None


def genome_differences(base, other, min_len):
    """Differences between two collinear genome versions as (b_start, b_end, o_start, o_end)."""
    al = mappy.Aligner(seq=base, preset="asm5", best_n=50)
    hits = [h for h in al.map(other) if h.strand == 1]
    if not hits:
        return [], 0.0
    hits.sort(key=lambda h: h.q_st)
    n = len(hits)
    best = [h.mlen for h in hits]
    prev = [-1] * n
    for i in range(n):
        for j in range(i):
            if hits[j].q_en - 1000 <= hits[i].q_st and hits[j].r_en - 1000 <= hits[i].r_st:
                if best[j] + hits[i].mlen > best[i]:
                    best[i] = best[j] + hits[i].mlen
                    prev[i] = j
    i = max(range(n), key=lambda k: best[k])
    chain = []
    while i >= 0:
        chain.append(hits[i])
        i = prev[i]
    chain.reverse()
    covered = sum(h.q_en - h.q_st for h in chain) / max(1, len(other))
    diffs = []
    for h1, h2 in zip(chain, chain[1:]):
        bs, be, os_, oe = h1.r_en, h2.r_st, h1.q_en, h2.q_st
        if be < bs or oe < os_:          # overlapping blocks: shrink both sides by the larger overlap
            ov = max(bs - be, os_ - oe, 0)
            bs, os_ = bs - ov, os_ - ov
        if max(be - bs, oe - os_) >= min_len:
            diffs.append((bs, be, os_, oe))
    for h in chain:
        rp, qp = h.r_st, h.q_st
        for ln, op in h.cigar:
            if op in (0, 7, 8):
                rp += ln
                qp += ln
            elif op == 1:
                if ln >= min_len:
                    diffs.append((rp, rp, qp, qp + ln))
                qp += ln
            elif op == 2:
                if ln >= min_len:
                    diffs.append((rp, rp + ln, qp, qp))
                rp += ln
            elif op == 3:
                rp += ln
    return sorted(diffs), covered


def classify_segments(seg_ref, seg_alt, left, right):
    if not seg_alt and seg_ref:
        return "DEL"
    if seg_alt and not seg_ref:
        if seg_alt == left[-len(seg_alt):] or seg_alt == right[:len(seg_alt)]:
            return "DUP"
        return "INS"
    if len(seg_ref) == len(seg_alt) and seg_alt == rc(seg_ref):
        return "INV"
    if len(seg_alt) == 2 * len(seg_ref) and seg_alt == seg_ref + seg_ref:
        return "DUP"
    return "COMPLEX"


def known_loci(paths, min_len, flank_base, anchor):
    genomes = []
    notes = []
    for p in paths:
        cs = read_fasta(p)
        if len(cs) != 1:
            notes.append(f"{p}: {len(cs)} sequences; only single-sequence genome versions are supported, skipped")
            continue
        genomes.append((os.path.splitext(os.path.basename(p))[0], cs[0]))
    if len(genomes) < 2:
        return [], notes + ["fewer than two usable genome versions"]
    base_name, base = genomes[0]
    brmap = RepeatMap([base])
    raw = []
    for vname, vc in genomes[1:]:
        seq, strand, rot = normalise_to(base.seq, vc.seq)
        if seq is None:
            notes.append(f"{vname}: could not be placed against {base_name}; skipped")
            continue
        diffs, cov = genome_differences(base.seq, seq, min_len)
        notes.append(f"{vname}: {len(diffs)} difference(s) >= {min_len} bp vs {base_name} "
                     f"(collinear chain covers {100 * cov:.2f}% of {vname}; orientation {strand}, rotation {rot})")
        for bs, be, os_, oe in diffs:
            raw.append((bs, be, vname, seq[os_:oe]))
    raw.sort()
    groups = []
    for r in raw:
        if groups and r[0] <= groups[-1]["be"] + 50:
            groups[-1]["items"].append(r)
            groups[-1]["be"] = max(groups[-1]["be"], r[1])
        else:
            groups.append({"bs": r[0], "be": r[1], "items": [r]})
    loci = []
    for g in groups:
        bs, be = g["bs"], g["be"]
        vsegs = []
        for it_bs, it_be, vname, seg in g["items"]:
            vsegs.append((vname, base.sub(bs, it_bs) + seg + base.sub(it_be, be)))
        lmax = max([be - bs] + [len(v) for _, v in vsegs])
        flank = flank_base + lmax // 2 + anchor
        left, right = base.sub(bs - flank, bs), base.sub(be, be + flank)
        if left is None or right is None:
            notes.append(f"difference at {base_name}:{bs + 1}-{be}: too close to a contig end for flanks; skipped")
            continue
        alleles = [Allele(base_name, left + base.sub(bs, be) + right, flank, flank + (be - bs), base_name)]
        seen = {alleles[0].seq: base_name}
        kinds = set()
        for vname, vseg in vsegs:
            s = left + vseg + right
            if s in seen:
                notes.append(f"{vname} carries the same sequence as {seen[s]} at {base_name}:{bs + 1}-{be}")
                continue
            seen[s] = vname
            alleles.append(Allele(vname, s, flank, flank + len(vseg), vname))
            kinds.add(classify_segments(base.sub(bs, be), vseg, left, right))
        if len(alleles) < 2:
            continue
        kind = kinds.pop() if len(kinds) == 1 else "COMPLEX"
        size = (be - bs) if kind == "INV" else max(abs(a.vlen - alleles[0].vlen) for a in alleles[1:])
        lid = f"known_{kind}{size}_{base_name}_{bs + 1}"
        desc = f"{base_name}:{bs + 1}-{be}" + "".join(f"; {a.name} segment {a.vlen} bp" for a in alleles)
        ml = brmap.mask_array(base.name, bs - flank, bs, base.len, base.circular)
        mr = brmap.mask_array(base.name, be, be + flank, base.len, base.circular)
        L = Locus(lid, "known", kind, alleles, desc, ml, mr)
        L.meta.update(base=base_name, base_start=bs, base_end=be)
        loci.append(L)
    return loci, notes


# ------------------------------------------------- Sniffles2 candidates
def run_sniffles(bam, ref, out_vcf, threads, min_len, af_min, logf):
    cmd = ["sniffles", "--input", bam, "--reference", ref, "--vcf", out_vcf, "--threads", str(threads),
           "--minsvlen", str(min_len), "--mosaic", "--mosaic-include-germline", "--mosaic-af-min", str(af_min),
           "--output-rnames", "--allow-overwrite"]
    log("[discovery] " + " ".join(cmd))
    with open(logf, "w") as fh:
        return subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)


def parse_vcf(path):
    recs = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = {}
            for kv in f[7].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info[k] = v
                else:
                    info[kv] = True
            recs.append({"chrom": f[0], "pos": int(f[1]), "id": f[2], "ref": f[3], "alt": f[4],
                         "filter": f[6], "info": info})
    return recs


def read_inserts(bam_path, chrom, pos, svlen, max_n=15, tol=300):
    """Inserted segments carried in single alignments (CIGAR I) near `pos`; deterministic subset by read name."""
    cands = []
    with pysam.AlignmentFile(bam_path) as fh:
        for a in fh.fetch(chrom, max(0, pos - tol), pos + tol):
            if a.is_secondary or a.is_supplementary or a.is_unmapped or a.query_sequence is None:
                continue
            rp_, qp = a.reference_start, 0
            for op, ln in a.cigartuples:
                if op == 1 and abs(rp_ - pos) <= tol and 0.8 * svlen <= ln <= 1.2 * svlen:
                    cands.append((name_hash(a.query_name), a.query_sequence[qp:qp + ln]))
                    break
                if op in (0, 7, 8):
                    rp_ += ln
                    qp += ln
                elif op in (1, 4):
                    qp += ln
                elif op in (2, 3):
                    rp_ += ln
    cands.sort()
    return [sq for _, sq in cands[:max_n]]


def seq_identity(a, b):
    """Best single-alignment matches of b on a, relative to the longer sequence."""
    hits = list(mappy.Aligner(seq=a, preset="map-ont").map(b))
    return max((h.mlen for h in hits), default=0) / max(len(a), len(b), 1)


def medoid(seqs):
    """The sequence most similar to all others (robust stand-in for a consensus of noisy reads)."""
    if len(seqs) <= 2:
        return seqs[0]
    score = [0.0] * len(seqs)
    for i, s in enumerate(seqs):
        al = mappy.Aligner(seq=s, preset="map-ont")
        for j, t in enumerate(seqs):
            if i != j:
                score[i] += max((h.mlen for h in al.map(t)), default=0) / max(len(s), len(t))
    return seqs[max(range(len(seqs)), key=lambda k: score[k])]


def discovery_loci(sample, vcf, contigs, rmap, min_len, flank_base, anchor, max_cand, bam_path=None):
    cmap = {c.name: c for c in contigs}
    cands = []
    skipped = {"BND/other": 0, "symbolic INS": 0, "origin-spanning": 0, "too small": 0}
    for r in parse_vcf(vcf):
        t = r["info"].get("SVTYPE", "")
        if t not in ("DEL", "INS", "INV", "DUP"):
            skipped["BND/other"] += 1
            continue
        c = cmap.get(r["chrom"])
        if c is None:
            continue
        try:
            svlen = abs(int(r["info"].get("SVLEN", "0")))
        except ValueError:
            svlen = 0
        if svlen < min_len:
            skipped["too small"] += 1
            continue
        start = r["pos"]
        end = start if t == "INS" else (start + svlen if t == "DEL" else int(r["info"].get("END", start + svlen)))
        if c.circular and (svlen >= 0.5 * c.len or (start < 5000 and end > c.len - 5000)):
            skipped["origin-spanning"] += 1
            continue
        alt_seq = None
        if t == "INS":
            alt = r["alt"]
            if alt.startswith("<"):
                skipped["symbolic INS"] += 1
                continue
            alt_seq = alt[1:] if alt[:1].upper() == r["ref"][:1].upper() else alt
        support = int(r["info"].get("SUPPORT", "0") or 0)
        af = float(r["info"].get("VAF", r["info"].get("AF", "0")) or 0)
        cands.append(dict(chrom=c.name, type=t, start=start, end=end, svlen=svlen, alt=alt_seq,
                          support=support, af=af, filter=r["filter"], id=r["id"]))
    cands.sort(key=lambda x: (x["chrom"], x["type"], x["start"]))
    merged = []
    for x in cands:
        if merged:
            y = merged[-1]
            if (y["chrom"] == x["chrom"] and y["type"] == x["type"] and abs(y["start"] - x["start"]) <= 100
                    and abs(y["svlen"] - x["svlen"]) <= 0.2 * max(y["svlen"], x["svlen"])):
                if x["support"] > y["support"]:
                    merged[-1] = x
                continue
        merged.append(x)
    merged.sort(key=lambda x: -x["support"])
    if len(merged) > max_cand:
        skipped["over candidate cap"] = len(merged) - max_cand
        merged = merged[:max_cand]
    loci = []
    for x in merged:
        c = cmap[x["chrom"]]
        flank = flank_base + x["svlen"] // 2 + anchor
        s, e = x["start"], x["end"]
        ins_note = ""
        if x["type"] == "DEL":
            alt = ""
        elif x["type"] == "INS":
            # Sniffles2 can report the right position and length with the wrong sequence (seen: a copy of
            # a region 19 kb away), so the insert is rebuilt from the reads that carry it in one alignment.
            alt = x["alt"]
            ins = read_inserts(bam_path, c.name, x["start"], x["svlen"]) if bam_path else []
            if len(ins) >= 3:
                m = medoid(ins)
                agree = seq_identity(m, x["alt"]) >= 0.9 if x["alt"] else False
                alt = m
                ins_note = (f"; insert = medoid of {len(ins)} read-carried inserts"
                            + ("" if agree else " (Sniffles2 sequence disagreed; not used)"))
            else:
                ins_note = "; insert sequence from Sniffles2 (fewer than 3 reads carry it in one alignment)"
        elif x["type"] == "INV":
            seg = c.sub(s, e)
            alt = rc(seg) if seg else None
        else:
            seg = c.sub(s, e)
            alt = seg + seg if seg else None
        if alt is None:
            continue
        lid = f"disc_{sample}_{x['type']}{x['svlen']}_{c.name}_{s + 1}"
        desc = (f"{sample} consensus {c.name}:{s + 1}-{max(e, s + 1)}; Sniffles2 {x['id']} SUPPORT={x['support']} "
                f"VAF={x['af']:.3f} FILTER={x['filter']}{ins_note}")
        drp = direct_repeat_pair(c, s, e) if x["type"] in ("DEL", "DUP") and e - s >= 300 else None
        if drp is not None:
            s2, e2 = drp[0], drp[2]
            if abs((s - s2) - (e - e2)) <= 10:      # breakpoints already one-copy-consistent
                drp = None
        if drp is None:
            L = make_locus_from_contig(c, rmap, lid, "discovered", x["type"], s, e, alt, flank, "consensus", x["type"], desc)
            if L is not None:
                L.meta["segments"] = [(L.alleles[0].vs, L.alleles[0].ve)]
        else:
            x0, x1 = min(s, s2) - flank, max(e, e2) + flank
            C = c.sub
            if x["type"] == "DEL":
                called = C(x0, s) + C(e, x1)
                hr = C(x0, s2) + C(e2, x1)
            else:
                called = C(x0, e) + C(s, e) + C(e, x1)
                hr = C(x0, e2) + C(s2, e2) + C(e2, x1)
            L = locus_from_full(c, rmap, lid, "discovered", x["type"], x0, x1,
                                [("REF", C(x0, x1), "consensus"), ("ALT", hr, f"{x['type']} one repeat copy at junction"),
                                 ("ALT_called", called, f"{x['type']} as called")], desc
                                + f"; breakpoints in direct repeat {c.name}:{drp[0] + 1}-{drp[1]} / {drp[2] + 1}-{drp[3]} "
                                  f"({drp[1] - drp[0]} bp): one-copy allele ALT added")
            if L is not None:
                L.meta["segments"] = [(s - x0, e - x0), (s2 - x0, e2 - x0)]
        if L is not None:
            L.meta["sniffles"] = x
            loci.append(L)
    return loci, skipped


# ------------------------------------------------- per-sample genotyping
class Sample:
    def __init__(self, name, consensus, bam, gff=None, fna=None):
        self.name, self.consensus, self.bam, self.gff, self.fna = name, consensus, bam, gff, fna


def locate(locus, aligner, cmap, probe=2000):
    """Place a locus in a sample consensus via unique hits of its outer flank ends."""
    a0 = locus.alleles[0].seq
    lp, rp = a0[:probe], a0[-probe:]
    hl = [h for h in aligner.map(lp) if h.q_en - h.q_st >= 0.9 * probe]
    hr = [h for h in aligner.map(rp) if h.q_en - h.q_st >= 0.9 * probe]
    if len(hl) != 1 or len(hr) != 1:
        return None, f"flank probes not unique (left {len(hl)} hits, right {len(hr)} hits)"
    hl, hr = hl[0], hr[0]
    if hl.ctg != hr.ctg or hl.strand != hr.strand:
        return None, "flank probes on different contigs or strands"
    c = cmap[hl.ctg]
    if hl.strand == 1:
        start, end = hl.r_st - hl.q_st, hr.r_en + (probe - hr.q_en)
    else:
        start, end = hr.r_st - (probe - hr.q_en), hl.r_en + hl.q_st
    span = end - start
    if span <= 0:
        if not c.circular:
            return None, "flank probes out of order on a linear contig"
        end += c.len
        span += c.len
    lens = [len(a.seq) for a in locus.alleles]
    if span > max(lens) * 1.2 + 1000 or span < min(lens) * 0.8 - 1000:
        return None, f"flank probes {span} bp apart, not consistent with allele lengths {lens}"
    match = [a.name for a in locus.alleles if abs(len(a.seq) - span) <= max(50, 0.02 * len(a.seq))]
    return {"contig": c.name, "start": start, "end": end, "strand": "+" if hl.strand == 1 else "-",
            "span": span, "consensus_allele": match[0] if len(match) == 1 else ""}, ""


def split_interval(contig, s, e):
    """Fetchable pieces of [s, e) on a (possibly circular) contig."""
    clen = contig.len
    if not contig.circular:
        return [(max(0, s), min(clen, e))]
    s2, e2 = s, e
    while s2 < 0:
        s2, e2 = s2 + clen, e2 + clen
    while s2 >= clen:
        s2, e2 = s2 - clen, e2 - clen
    if e2 <= clen:
        return [(s2, e2)]
    return [(s2, clen), (0, e2 - clen)]


def anchor_regions(L, loc):
    """Consensus coordinates of the locus anchors (shared by every allele), or None."""
    a = next((x for x in L.alleles if x.name == loc.get("consensus_allele")), None)
    if a is None or L.A_left is None:
        return None
    if loc["strand"] == "+":
        bl, br = loc["start"] + a.vs, loc["start"] + a.ve
        return (bl - L.A_left, bl), (br, br + L.A_right)
    bl, br = loc["end"] - a.vs, loc["end"] - a.ve
    return (bl, bl + L.A_left), (br - L.A_right, br)


def fetch_reads(bam, contig, start, end, min_len, max_reads, regions=None, slack=50):
    """Reads that can span the locus.

    With `regions` (the two anchor intervals in consensus coordinates), a read
    is kept if its primary or supplementary alignments cover both anchors;
    every allele shares the anchors, so any read able to span any allele's
    window passes this test whichever allele it carries. Without regions,
    all primary alignments overlapping [start, end) are kept. The read
    sequence always comes from the primary record (it keeps soft clips)."""
    prim = {}
    stats = {"primary_elsewhere": 0}
    if regions is None:
        for s, e in split_interval(contig, start, end):
            if e <= s:
                continue
            for r in bam.fetch(contig.name, s, e):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                if r.query_name not in prim:
                    prim[r.query_name] = (r.query_sequence, r.is_reverse)
        names = list(prim)
    else:
        need, got = {}, {}
        for ri, (s, e) in enumerate(regions):
            pieces = split_interval(contig, s, e)
            need[ri] = len(pieces)
            for pi, (fs, fe) in enumerate(pieces):
                if fe <= fs:
                    continue
                for r in bam.fetch(contig.name, fs, fe):
                    if r.is_unmapped or r.is_secondary:
                        continue
                    nm = r.query_name
                    if not r.is_supplementary and nm not in prim:
                        prim[nm] = (r.query_sequence, r.is_reverse)
                    if r.reference_start <= fs + slack and r.reference_end >= fe - slack:
                        got.setdefault(nm, set()).add((ri, pi))
        names = []
        for nm, cov in got.items():
            if all(sum(1 for (ri, _) in cov if ri == k) == need[k] for k in need):
                if nm in prim:
                    names.append(nm)
                else:
                    stats["primary_elsewhere"] += 1
    reads = [(nm, prim[nm][0], prim[nm][1]) for nm in names
             if prim[nm][0] is not None and len(prim[nm][0]) >= min_len]
    n_elig = len(reads)
    if len(reads) > max_reads:
        reads.sort(key=lambda x: name_hash(x[0]))
        reads = reads[:max_reads]
    return reads, n_elig, stats


_tls = threading.local()


def span_check(hit, win, max_indel, min_ident):
    """Identity of `hit` if it covers `win` with no internal indel >= max_indel, else None."""
    if hit.r_st > win[0] or hit.r_en < win[1]:
        return None
    if hit.mlen < min_ident * hit.blen:
        return None
    pos = hit.r_st
    for ln, op in hit.cigar:
        if op in (0, 7, 8, 3):
            pos += ln
        elif op == 2:
            if ln >= max_indel and pos < win[1] and pos + ln > win[0]:
                return None
            pos += ln
        elif op == 1:
            if ln >= max_indel and win[0] <= pos <= win[1]:
                return None
    return hit.mlen / hit.blen


def diverges(h, alen, rlen, min_clip=500, tol=200):
    """True if the read leaves the allele sequence inside it: >= min_clip unaligned read bases on a side
    where the allele continues for > tol bp. Such reads (construct-borne junctions such as homology-arm
    plasmids, chimeric reads, reads carrying another nearby difference) are not counted for any allele."""
    before, after = h.q_st, rlen - h.q_en
    if h.strand == -1:
        before, after = after, before
    return (h.r_st > tol and before >= min_clip) or (h.r_en < alen - tol and after >= min_clip)


def genotype(locus, reads, aligners, min_ident=0.85, margin=0.01, threads=4):
    wins = [a.win for a in locus.alleles]

    def one(rd):
        nm, seq, rev = rd
        buf = getattr(_tls, "buf", None)
        if buf is None:
            buf = _tls.buf = mappy.ThreadBuffer()
        scores, div = [], False
        for k, al in enumerate(aligners):
            best = None
            for h in al.map(seq, buf=buf):
                s = span_check(h, wins[k], locus.max_indel, min_ident)
                if s is None:
                    continue
                if diverges(h, len(locus.alleles[k].seq), len(seq)):
                    div = True
                    continue
                if best is None or s > best:
                    best = s
            scores.append(best)
        ok = sorted(((s, k) for k, s in enumerate(scores) if s is not None), reverse=True)
        if not ok:
            return nm, rev, "divergent" if div else None
        if len(ok) >= 2 and ok[0][0] - ok[1][0] < margin:
            return nm, rev, "ambiguous"
        return nm, rev, ok[0][1]

    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        res = list(ex.map(one, reads))
    n = len(locus.alleles)
    counts, fwd, names = [0] * n, [0] * n, [[] for _ in range(n)]
    amb = none = divg = 0
    for nm, rev, k in res:
        if k is None:
            none += 1
        elif k == "divergent":
            divg += 1
        elif k == "ambiguous":
            amb += 1
        else:
            counts[k] += 1
            fwd[k] += 0 if rev else 1
            names[k].append(nm)
    return {"counts": counts, "fwd": fwd, "names": names, "ambiguous": amb, "unspanned": none, "divergent": divg}


def coverage_ratio(bam, contig, seg_start, seg_end, rmap, n=25, gap=2000, width=10000):
    """Median depth inside [seg_start, seg_end) over median depth at unique positions 2-12 kb outside."""
    def depth(p):
        p %= contig.len
        return sum(1 for r in bam.fetch(contig.name, p, p + 1) if not (r.is_unmapped or r.is_secondary))
    L = seg_end - seg_start
    if L < 1000:
        return None
    inside = [seg_start + int((i + 0.5) * L / n) for i in range(n)]
    outside = []
    for i in range(n):
        off = gap + int((i + 0.5) * width / n)
        for p in (seg_start - off, seg_end + off):
            if contig.circular or 0 <= p < contig.len:
                q = p % contig.len
                if rmap.masked_fraction(contig.name, q, q + 1) == 0:
                    outside.append(p)
    if len(outside) < 5:
        return None
    din = sorted(depth(p) for p in inside)
    dout = sorted(depth(p) for p in outside)
    mo = dout[len(dout) // 2]
    return round(din[len(din) // 2] / mo, 3) if mo else None


# ------------------------------------------------- breakpoint-junction genotyping
# For events whose full allele window is longer than the reads, each read is
# instead assigned to one of four equal-geometry junction alleles. Segment
# events (DEL / INS / DUP / excision): PRESENT_LEFT (left flank|segment start),
# PRESENT_RIGHT (segment end|right flank), END_TO_START (segment end|segment
# start: tandem copy or circular form) and ABSENT (left flank|right flank).
# Inversions: REF_LEFT, INV_LEFT, REF_RIGHT, INV_RIGHT. Breakpoints are taken
# as reported (read-supported), so repeats adjacent to a junction are kept
# exactly as the reads show them; windows extend past such repeats until a
# unique anchor is reached on both sides, and every allele uses the same window
# geometry, so no allele is favoured by length.
SEG_JUNCTIONS = ("PRESENT_LEFT", "PRESENT_RIGHT", "END_TO_START", "ABSENT")
INV_JUNCTIONS = ("REF_LEFT", "INV_LEFT", "REF_RIGHT", "INV_RIGHT")


def piece_mask(piece, cal, chunk=1000, step=250, min_frac=0.9, min_ident=0.85):
    """Positions of `piece` inside sequence present >= 2 times in the sample consensus."""
    m = [False] * len(piece)
    if len(piece) < chunk:
        return m
    starts = list(range(0, len(piece) - chunk + 1, step))
    if starts[-1] != len(piece) - chunk:
        starts.append(len(piece) - chunk)
    for st in starts:
        q = piece[st:st + chunk]
        hits = [h for h in cal.map(q) if h.q_en - h.q_st >= min_frac * chunk and h.mlen >= min_ident * h.blen]
        if len(hits) >= 2:
            for p in range(st, st + chunk):
                m[p] = True
    return m


def unique_reach(mask, left_of_junction, anchor):
    """bp from the junction needed to include `anchor` consecutive unique bp of the piece, or None."""
    n, run = len(mask), 0
    for d in range(1, n + 1):
        run = 0 if mask[n - d if left_of_junction else d - 1] else run + 1
        if run >= anchor:
            return d
    return None


def unique_hit(part, cal, min_frac=0.8):
    hits = [h for h in cal.map(part) if h.q_en - h.q_st >= min_frac * len(part)]
    if len(hits) != 1:
        return None
    return (hits[0].ctg, hits[0].r_st, hits[0].r_en)


def junction_context(L):
    """How to build junction alleles for a locus: dict(model, S, s, e, role) or (None, reason).

    role tells how the junction estimates map onto the locus alleles:
      seg_ref   segment = alleles[0] variable part (DEL/DUP/INV-like; consensus or base carries it)
      seg_alt   segment = alleles[1] variable part (insertion carried by the second allele)
      seg_X     segment carried by allele index X (known presence/absence pairs)"""
    a = L.alleles
    if L.kind == "INV":
        return {"model": "INV", "S": a[0].seq, "s": a[0].vs, "e": a[0].ve, "seg_allele": 0}, ""
    v0, v1 = a[0].vlen, a[1].vlen if len(a) > 1 else 0
    if L.kind == "DUP" and L.source == "known":
        if len(a) != 2:
            return None, "junction mode supports two-allele loci only"
        base, dup = (a[0], a[1]) if v0 < v1 else (a[1], a[0])
        seg = dup.seq[dup.vs:dup.ve]
        if len(seg) == 2 * base.vlen and seg == base.seq[base.vs:base.ve] * 2:
            s, e = base.vs, base.ve
        elif base.vlen == 0 and seg == base.seq[base.ve:base.ve + len(seg)]:
            s, e = base.vs, base.vs + len(seg)
        elif base.vlen == 0 and seg == base.seq[base.vs - len(seg):base.vs]:
            s, e = base.vs - len(seg), base.vs
        else:
            return None, "duplication not representable as a tandem copy"
        return {"model": "SEG", "S": base.seq, "s": s, "e": e, "seg_allele": a.index(base), "dup_allele": a.index(dup)}, ""
    segs = L.meta.get("segments")
    if L.kind in ("DUP", "DEL") and L.source != "known" and segs:
        ctx = {"model": "SEG", "S": a[0].seq, "s": segs[0][0], "e": segs[0][1], "seg_allele": 0, "alt_segments": segs[1:]}
        if L.kind == "DUP":
            ctx["dup_allele"] = 1
        return ctx, ""
    if L.kind == "DUP":          # null loci: alleles[0] = consensus with one copy of [vs, ve)
        return {"model": "SEG", "S": a[0].seq, "s": a[0].vs, "e": a[0].ve, "seg_allele": 0, "dup_allele": 1}, ""
    if len(a) != 2:
        return None, "junction mode supports two-allele loci only"
    if min(v0, v1) > 50:
        return None, "both alleles carry long sequence; not a presence/absence segment"
    k = 0 if v0 >= v1 else 1
    return {"model": "SEG", "S": a[k].seq, "s": a[k].vs, "e": a[k].ve, "seg_allele": k}, ""


def junction_qpos(h, rpos):
    """Read coordinate aligned to allele position `rpos` in mappy hit `h` (original read orientation)."""
    r, q = h.r_st, 0
    for ln, op in h.cigar:
        if op in (0, 7, 8):
            if r + ln >= rpos:
                q += rpos - r
                break
            r += ln
            q += ln
        elif op in (2, 3):
            if r + ln >= rpos:
                break
            r += ln
        elif op == 1:
            q += ln
    return h.q_st + q if h.strand == 1 else h.q_en - q


def build_junctions(L, cal, flank, anchor):
    """Junction alleles (equal window geometry) for locus L, or (None, reason)."""
    ctx, why = junction_context(L)
    if ctx is None:
        return None, why
    S, s, e = ctx["S"], ctx["s"], ctx["e"]
    segs = [(s, e)] + [tuple(x) for x in ctx.get("alt_segments", [])]
    F = min([flank] + [v for s_, e_ in segs for v in (s_, len(S) - e_)])
    Fs = min([flank] + [e_ - s_ for s_, e_ in segs])
    if F < anchor + 200 or Fs < anchor + 50:
        return None, f"flank {F} bp / segment {e - s} bp too short for {anchor}-bp junction anchors"
    P = {}
    for i, (s_, e_) in enumerate(segs):
        x = "" if i == 0 else str(i)
        P.update({"L" + x: S[s_ - F:s_], "Ms" + x: S[s_:s_ + Fs], "Me" + x: S[e_ - Fs:e_], "R" + x: S[e_:e_ + F]})
    M = {k: piece_mask(v, cal) for k, v in P.items()}
    P["rcMe"], M["rcMe"] = rc(P["Me"]), M["Me"][::-1]
    P["rcMs"], M["rcMs"] = rc(P["Ms"]), M["Ms"][::-1]
    if ctx["model"] == "SEG":
        defs = [("PRESENT_LEFT", "L", "Ms"), ("PRESENT_RIGHT", "Me", "R"), ("END_TO_START", "Me", "Ms"), ("ABSENT", "L", "R")]
        # alternative breakpoint hypotheses (e.g. one repeat copy at the junction): same junction types
        for i in range(1, len(segs)):
            defs += [(f"END_TO_START:alt{i}", f"Me{i}", f"Ms{i}"), (f"ABSENT:alt{i}", f"L{i}", f"R{i}")]
    else:
        defs = [("REF_LEFT", "L", "Ms"), ("INV_LEFT", "L", "rcMe"), ("REF_RIGHT", "Me", "R"), ("INV_RIGHT", "rcMs", "R")]
    reach = []
    for name, lk, rk in defs:
        rl, rr = unique_reach(M[lk], True, anchor), unique_reach(M[rk], False, anchor)
        if rl is None or rr is None:
            return None, f"no {anchor}-bp unique anchor within {flank} bp of the {name} junction"
        reach.append((rl, rr))
    RL, RR = max(r[0] for r in reach), max(r[1] for r in reach)
    if any(len(P[lk]) < RL + 50 or len(P[rk]) < RR + 50 for _, lk, rk in defs):
        return None, "junction pieces shorter than the repeat-extended anchors"
    alleles = []
    for name, lk, rk in defs:
        al = Allele(name, P[lk] + P[rk], len(P[lk]), len(P[lk]))
        al.win = (len(P[lk]) - RL, len(P[lk]) + RR)
        alleles.append(al)
    # anchor intervals in the sample consensus used to fetch reads (absent pieces are simply not found)
    regions, found = [], {}
    for k in [k for k in P if not k.startswith("rc")]:
        part = P[k][len(P[k]) - RL:] if k[0] in ("L",) or k.startswith("Me") else P[k][:RR]
        h = unique_hit(part, cal)
        found[k] = h is not None
        if h is not None and h not in regions:
            regions.append(h)
    return {"model": ctx["model"], "ctx": ctx, "alleles": alleles, "RL": RL, "RR": RR, "regions": regions,
            "found": found, "W": RL + RR}, ""


def fetch_any(bam, regions, min_len, max_reads, slack=50):
    """Reads with a primary or supplementary alignment covering any of the regions."""
    prim, got = {}, set()
    for ctg, s, e in regions:
        for r in bam.fetch(ctg, s, e):
            if r.is_unmapped or r.is_secondary:
                continue
            nm = r.query_name
            if not r.is_supplementary and nm not in prim:
                prim[nm] = (r.query_sequence, r.is_reverse)
            if r.reference_start <= s + slack and r.reference_end >= e - slack:
                got.add(nm)
    names = [nm for nm in got if nm in prim]
    elsewhere = len(got) - len(names)
    reads = [(nm, prim[nm][0], prim[nm][1]) for nm in names if prim[nm][0] is not None and len(prim[nm][0]) >= min_len]
    n_el = len(reads)
    if len(reads) > max_reads:
        reads.sort(key=lambda x: name_hash(x[0]))
        reads = reads[:max_reads]
    return reads, n_el, elsewhere


def genotype_junctions(J, reads, aligners, min_ident=0.85, margin=0.01, threads=4, max_indel=200):
    """Per read region: each region of a read that spans a junction window supports one junction allele.
    Different regions of the same read may support different junctions (e.g. both ends of a short segment)."""
    alleles = J["alleles"]
    W = J["W"]

    def one(rd):
        nm, seq, rev = rd
        buf = getattr(_tls, "buf", None)
        if buf is None:
            buf = _tls.buf = mappy.ThreadBuffer()
        spans, ndiv = [], 0
        for k, al in enumerate(aligners):
            for h in al.map(seq, buf=buf):
                sc = span_check(h, alleles[k].win, max_indel, min_ident)
                if sc is None:
                    continue
                if diverges(h, len(alleles[k].seq), len(seq)):
                    ndiv += 1
                    continue
                # cluster on the read coordinate of the junction itself, not on the whole hit
                spans.append((sc, k, junction_qpos(h, alleles[k].vs)))
        spans.sort(reverse=True)
        clusters = []          # [best_ident, best_k, junction_read_pos, ambiguous]
        for sc, k, qp in spans:
            for cl in clusters:
                if abs(qp - cl[2]) < 0.5 * W:
                    if k != cl[1] and cl[0] - sc < margin:
                        cl[3] = True
                    break
            else:
                clusters.append([sc, k, qp, False])
        return nm, rev, [(c[1], c[3]) for c in clusters], ndiv

    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        res = list(ex.map(one, reads))
    n = len(alleles)
    counts, fwd, names = [0] * n, [0] * n, [[] for _ in range(n)]
    amb = none = divg = 0
    for nm, rev, cls, ndiv in res:
        if not cls:
            if ndiv:
                divg += 1
            else:
                none += 1
        for k, is_amb in cls:
            if is_amb:
                amb += 1
                continue
            counts[k] += 1
            fwd[k] += 0 if rev else 1
            names[k].append(nm)
    # per junction type (alternative breakpoint hypotheses of one type are summed)
    types = []
    for a in alleles:
        t = a.name.split(":")[0]
        if t not in types:
            types.append(t)
    tcounts = [sum(c for a, c in zip(alleles, counts) if a.name.split(":")[0] == t) for t in types]
    return {"counts": counts, "fwd": fwd, "names": names, "ambiguous": amb, "unspanned": none, "divergent": divg,
            "types": types, "type_counts": tcounts}


def junction_estimates(model, counts, seg_in_consensus):
    """Point estimates and Wilson 95% intervals from junction counts.

    SEG: segment_absent = ABSENT / (ABSENT + mean(PRESENT_LEFT, PRESENT_RIGHT)) (fraction of genomes
    lacking the segment); end_to_start_per_genome = END_TO_START / (same denominator), reported only
    when the consensus carries the segment at this locus (otherwise those reads cannot be recruited).
    INV: inverted = inverted-junction observations / all junction observations."""
    out = {}
    if model == "SEG":
        pl, pr, es, ab = counts
        ref = (pl + pr) / 2.0
        n = ab + int(round(ref))
        x = ab / (ab + ref) if ab + ref > 0 else None
        lo, hi = wilson(ab, n)
        out["segment_absent"] = (x, lo, hi)
        out["segment_present"] = (None if x is None else 1 - x, 1 - hi, 1 - lo)
        if seg_in_consensus:
            n2 = es + n
            lo2, hi2 = wilson(es, n2)
            conv = lambda p: p / (1 - p) if p < 1 else float("inf")
            out["end_to_start_per_genome"] = (es / (ab + ref) if ab + ref > 0 else None, conv(lo2), conv(hi2))
        out["n"] = n
    else:
        rl, il, rr, ir = counts
        n = rl + il + rr + ir
        v = (il + ir) / n if n else None
        lo, hi = wilson(il + ir, n)
        out["inverted"] = (v, lo, hi)
        out["n"] = n
    return out


def null_loci(sample, contigs, rmap, signature, nloci, flank_base, anchor, avoid, seed):
    """Size-matched fake events at random unique positions of the consensus."""
    kind, size = signature
    rng = random.Random(f"{seed}:{sample}:{kind}:{size}")
    flank = flank_base + size // 2 + anchor
    eligible = [c for c in contigs if c.len > 2 * (size + 2 * flank) + 40000]
    if not eligible:
        return []
    weights = [c.len for c in eligible]
    out, tries = [], 0
    while len(out) < nloci and tries < nloci * 50:
        tries += 1
        c = rng.choices(eligible, weights)[0]
        lo, hi = flank + 10000, c.len - size - flank - 10000
        if hi <= lo:
            continue
        x = rng.randrange(lo, hi)
        s, e = x, x + (size if kind in ("DEL", "INV", "DUP") else 0)
        if rmap.masked_fraction(c.name, s - flank, e + flank) > 0.05:
            continue
        if any(a_c == c.name and s - 30000 < a_e and e + 30000 > a_s for a_c, a_s, a_e in avoid):
            continue
        if kind == "DEL":
            alt = ""
        elif kind == "INV":
            alt = rc(c.sub(s, e))
        elif kind == "DUP":
            seg = c.sub(s, e)
            alt = seg + seg
        else:
            d = rng.randrange(flank, c.len - size - flank)
            if abs(d - x) < 50000 or rmap.masked_fraction(c.name, d, d + size) > 0.05:
                continue
            alt = c.sub(d, d + size)
        L = make_locus_from_contig(c, rmap, f"null_{kind}{size}_{len(out) + 1}", "null", kind, s, e, alt, flank,
                                   "consensus", "fake", f"{c.name}:{s + 1}")
        if L is not None:
            out.append(L)
    return out


def contig_name_map(cons_contigs, fna_path):
    """Consensus contig name -> annotation contig name, matched by sequence."""
    if not fna_path or not os.path.exists(fna_path):
        return {}
    by_seq = {hashlib.md5(c.seq.encode()).hexdigest(): c.name for c in read_fasta(fna_path)}
    return {c.name: by_seq[h] for c in cons_contigs if (h := hashlib.md5(c.seq.encode()).hexdigest()) in by_seq}


def load_gff(path):
    feats = {}
    if not path or not os.path.exists(path):
        return feats
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] in ("gene", "region", "source"):
                continue
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            prod = attr.get("product", "").replace("%2C", ",").replace("%3B", ";")
            feats.setdefault(f[0], []).append((int(f[3]), int(f[4]), f[2], attr.get("locus_tag", attr.get("ID", "")),
                                               attr.get("gene", ""), prod))
    return feats


def feat_label(f):
    """locus_tag(gene):product, the form used in every gene column."""
    return f"{f[3]}{'(' + f[4] + ')' if f[4] else ''}:{f[5]}"


def genes_in(feats, seqid, start1, end1):
    """Features overlapping start1..end1 (1-based, inclusive). start1 > end1 is a segment
    that runs across the origin of a circular contig (start1..contig end, then 1..end1)."""
    out = []
    wrap = start1 > end1
    for f in feats.get(seqid, []):
        s, e = f[0], f[1]
        if (e >= start1 or s <= end1) if wrap else (e >= start1 and s <= end1):
            out.append(feat_label(f))
    return out


def _covers(f, pos, clen):
    """Feature f contains 1-based pos; a feature with end > contig length crosses the origin."""
    return f[0] <= pos <= f[1] or (f[1] > clen and pos <= f[1] - clen)


def _cut(f, pos, clen, side):
    """Feature f (which contains pos) extends beyond a segment that starts or ends at pos."""
    if f[0] <= pos <= f[1]:
        return f[0] < pos if side == "start" else f[1] > pos
    return True if side == "start" else (f[1] - clen) > pos


def _flanks(fl, lo, hi, clen, circular):
    """Nearest feature ending at or before lo and nearest starting at or after hi, with the
    number of bases between it and the site (lo..hi exclusive); wraps on circular contigs."""
    left = [f for f in fl if f[1] <= lo]
    right = [f for f in fl if f[0] >= hi]
    fL = fR = None
    if left:
        fL = max(left, key=lambda f: (f[1], f[0], f[3]))
        gL = lo - fL[1]
    elif circular and fl:
        fL = max(fl, key=lambda f: (f[1], f[0], f[3]))
        gL = lo + (clen - fL[1])
    if right:
        fR = min(right, key=lambda f: (f[0], f[1], f[3]))
        gR = fR[0] - hi
    elif circular and fl:
        fR = min(fl, key=lambda f: (f[0], f[1], f[3]))
        gR = (clen - hi) + fR[0]
    lt = f"left {feat_label(fL)}, {gL} bp" if fL is not None else "left none"
    rt = f"right {feat_label(fR)}, {gR} bp" if fR is not None else "right none"
    return f"{lt}; {rt}"


def segment_genes(feats, seqid, ctg, S, E, clen, circular):
    """Genes at both ends of a variable segment S..E present in the consensus (1-based,
    inclusive; S > E runs across the origin). [cut]: the gene extends beyond the segment."""
    fl = feats.get(seqid, [])

    def at(pos, side):
        inside = sorted((f for f in fl if _covers(f, pos, clen)), key=lambda f: (f[0], f[1], f[3]))
        if inside:
            return "in " + ", ".join(feat_label(f) + (" [cut]" if _cut(f, pos, clen, side) else "") for f in inside)
        return "intergenic (" + _flanks(fl, pos - 1, pos + 1, clen, circular) + ")"
    return f"start {ctg}:{S} {at(S, 'start')} | end {ctg}:{E} {at(E, 'end')}"


def junction_genes(feats, seqid, ctg, p, q, clen, circular):
    """Genes at a junction between consensus bases p and q (q = p + 1, or 1 after the last
    base of a circular contig): where a sequence absent from the consensus joins it."""
    fl = feats.get(seqid, [])
    inside = sorted((f for f in fl if _covers(f, p, clen) and _covers(f, q, clen)), key=lambda f: (f[0], f[1], f[3]))
    if inside:
        ctx = "in " + ", ".join(feat_label(f) for f in inside)
    else:
        ctx = "intergenic (" + _flanks(fl, p, q, clen, circular) + ")"
    return f"junction {ctg}:{p}^{q} {ctx}"


# ------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True,
                    help="TSV: sample, consensus FASTA, sorted+indexed BAM, [annotation GFF3], [annotation FASTA]")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--variant-genomes", nargs="*", default=[],
                    help="FASTA files, one genome version each; the first is the base for comparisons")
    ap.add_argument("--discovery", action="store_true", help="also call candidates de novo with Sniffles2")
    ap.add_argument("--min-len", type=int, default=500)
    ap.add_argument("--anchor", type=int, default=1000)
    ap.add_argument("--flank-extra", type=int, default=8000, help="flank beyond anchor + half the event length")
    ap.add_argument("--min-alt-reads", type=int, default=5)
    ap.add_argument("--min-fraction", type=float, default=0.01)
    ap.add_argument("--min-spanning", type=int, default=20)
    ap.add_argument("--null-loci", type=int, default=20)
    ap.add_argument("--max-reads", type=int, default=3000)
    ap.add_argument("--max-candidates", type=int, default=200)
    ap.add_argument("--min-identity", type=float, default=0.85)
    ap.add_argument("--junction-mode", choices=("auto", "always", "never"), default="auto",
                    help="breakpoint-junction genotyping: auto = for events too long for spanning reads")
    ap.add_argument("--max-span-window", type=int, default=20000,
                    help="loci whose allele windows are longer than this are genotyped primarily at their breakpoint "
                         "junctions (few reads span such windows); spanning-read results are still reported")
    ap.add_argument("--junction-flank", type=int, default=8000,
                    help="bp of sequence on each side of a junction available for repeat-extended anchors")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    t0 = time.time()
    log(f"sv_genotype {VERSION}: {' '.join(sys.argv[1:])}")

    samples = []
    with open(args.samples) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t") + ["", ""]
            samples.append(Sample(f[0], f[1], f[2], f[3] or None, f[4] or None))

    known, notes = (known_loci(args.variant_genomes, args.min_len, args.flank_extra, args.anchor)
                    if args.variant_genomes else ([], []))
    for n in notes:
        log("[known] " + n)

    sdata = {}
    all_loci = list(known)
    for s in samples:
        sd = {"contigs": read_fasta(s.consensus)}
        sd["cmap"] = {c.name: c for c in sd["contigs"]}
        log(f"[{s.name}] repeat map of {len(sd['contigs'])} contig(s)")
        sd["rmap"] = RepeatMap(sd["contigs"])
        sd["outdir"] = os.path.join(args.outdir, s.name)
        os.makedirs(sd["outdir"], exist_ok=True)
        sd["disc_note"] = ""
        if args.discovery:
            ref = os.path.join(sd["outdir"], "consensus.fasta")
            with open(ref, "w") as fh:
                for c in sd["contigs"]:
                    fh.write(f">{c.name} {c.comment}\n{c.seq}\n")
            pysam.faidx(ref)
            vcf = os.path.join(sd["outdir"], "sniffles.vcf")
            rcode = run_sniffles(s.bam, ref, vcf, args.threads, args.min_len, args.min_fraction,
                                 os.path.join(sd["outdir"], "sniffles.log"))
            if rcode != 0 or not os.path.exists(vcf):
                log(f"[{s.name}] Sniffles2 failed (rc={rcode}); see sniffles.log")
                sd["disc_note"] = f"Sniffles2 rc={rcode}"
            else:
                loci, skipped = discovery_loci(s.name, vcf, sd["contigs"], sd["rmap"], args.min_len,
                                               args.flank_extra, args.anchor, args.max_candidates, bam_path=s.bam)
                log(f"[{s.name}] Sniffles2: {len(loci)} candidate loci; skipped {skipped}")
                sd["disc_note"] = f"{len(loci)} candidates; skipped {skipped}"
                all_loci.extend(loci)
        sdata[s.name] = sd

    for L in all_loci:
        compute_windows(L, args.anchor)

    rows, null_rows = [], []
    summary = {"version": VERSION, "args": vars(args), "known_notes": notes, "samples": {}}
    for s in samples:
        sd = sdata[s.name]
        bam = pysam.AlignmentFile(s.bam)
        cal = mappy.Aligner(s.consensus, preset="asm10", best_n=5)
        name_map = contig_name_map(sd["contigs"], s.fna)
        feats = load_gff(s.gff)
        located, merged = [], []
        for L in all_loci:
            loc, why = locate(L, cal, sd["cmap"])
            if loc is None:
                if L.source == "known" or L.id.startswith(f"disc_{s.name}_"):
                    log(f"[{s.name}] {L.id}: not placed ({why})")
                continue
            same = next((o for o, o_loc in located if o_loc["contig"] == loc["contig"]
                         and abs(o_loc["start"] - loc["start"]) < 200 and o.kind == L.kind
                         and abs(o.lmax - L.lmax) <= 0.2 * max(o.lmax, L.lmax, 1)), None)
            if same is not None:
                merged.append([L.id, same.id])
                if L.source == "discovered" and same.source == "known":
                    log(f"[{s.name}] {L.id}: same event as {same.id}; genotyped once, as the known locus")
                continue
            located.append((L, loc))
        log(f"[{s.name}] {len(located)} locus/loci placed in the consensus")

        def sig(L, loc):
            k = L.kind
            if k in ("DEL", "INS", "COMPLEX"):
                ca = loc.get("consensus_allele", "")
                cons_len = next((a.vlen for a in L.alleles if a.name == ca), L.lmax)
                k = "DEL" if cons_len == L.lmax else "INS"
            size = L.lmax if L.kind == "INV" else L.dlen
            # null events use fixed size bins (upper bin edge), so the null work stays bounded
            for b in NULL_BINS:
                if size <= b:
                    return (k, b)
            return (k, NULL_BINS[-1])

        avoid = [(loc["contig"], loc["start"], loc["end"]) for _, loc in located]
        floors = {}
        for sg in sorted({sig(L, loc) for L, loc in located if L.W is not None and L.W <= 60000}):
            fr = []
            for NL in null_loci(s.name, sd["contigs"], sd["rmap"], sg, args.null_loci, args.flank_extra, args.anchor,
                                avoid, args.seed):
                if not compute_windows(NL, args.anchor):
                    continue
                c = sd["cmap"][NL.meta["contig"]]
                st, en = NL.meta["start"], NL.meta["end"]
                regs = ((st - NL.A_left, st), (en, en + NL.A_right))
                reads, n_el, _ = fetch_reads(bam, c, st - NL.FL, en + NL.FR, NL.W, args.max_reads, regions=regs)
                al = [mappy.Aligner(seq=a.seq, preset="map-ont", best_n=5) for a in NL.alleles]
                g = genotype(NL, reads, al, args.min_identity, threads=args.threads)
                n = sum(g["counts"])
                f_alt = g["counts"][1] / n if n else None
                null_rows.append([s.name, f"{sg[0]}{sg[1]}", NL.id, NL.desc, NL.W, n_el, len(reads), g["counts"][0],
                                  g["counts"][1], g["ambiguous"], g["unspanned"], "" if f_alt is None else f"{f_alt:.5f}"])
                if n >= args.min_spanning:
                    fr.append(f_alt)
            floors[sg] = max([args.min_fraction] + fr)
            log(f"[{s.name}] null {sg[0]}{sg[1]}: {len(fr)} usable loci, max fake-allele fraction "
                f"{(max(fr) if fr else float('nan')):.4f} -> floor {floors[sg]:.4f}")

        jfloors = {}

        def junction_floor(jsg):
            """Max event-indicating junction fraction over size-matched null loci (cached per signature)."""
            if jsg in jfloors:
                return jfloors[jsg]
            kind, b = jsg
            fr = []
            for NL in null_loci(s.name, sd["contigs"], sd["rmap"], (kind, b), args.null_loci, args.flank_extra,
                                args.anchor, avoid, args.seed):
                J, _ = build_junctions(NL, cal, args.junction_flank, args.anchor)
                if J is None:
                    continue
                jr, n_el, _ = fetch_any(bam, J["regions"], J["W"], args.max_reads)
                jal = [mappy.Aligner(seq=a.seq, preset="map-ont", best_n=5) for a in J["alleles"]]
                g = genotype_junctions(J, jr, jal, args.min_identity, threads=args.threads)
                cnt = g["type_counts"]
                est = junction_estimates(J["model"], cnt, kind != "INS")
                if J["model"] == "INV":
                    stat, ev, rf = est["inverted"][0], cnt[1] + cnt[3], cnt[0] + cnt[2]
                elif kind == "INS":
                    stat, ev, rf = est["segment_present"][0], cnt[0] + cnt[1], cnt[3]
                else:
                    vals = [est["segment_absent"][0], est.get("end_to_start_per_genome", (None,))[0]]
                    vals = [v for v in vals if v is not None]
                    stat, ev, rf = (max(vals) if vals else None), cnt[2] + cnt[3], cnt[0] + cnt[1]
                null_rows.append([s.name, f"junction_{kind}{b}", NL.id, NL.desc, J["W"], n_el, len(jr), rf, ev,
                                  g["ambiguous"], g["unspanned"], "" if stat is None else f"{stat:.5f}"])
                if est["n"] >= args.min_spanning and stat is not None:
                    fr.append(stat)
            jfloors[jsg] = max([args.min_fraction] + fr)
            log(f"[{s.name}] junction null {kind}{b}: {len(fr)} usable loci, max fake fraction "
                f"{(max(fr) if fr else float('nan')):.4f} -> floor {jfloors[jsg]:.4f}")
            return jfloors[jsg]

        def run_junction(L, loc):
            J, why = build_junctions(L, cal, args.junction_flank, args.anchor)
            if J is None:
                return {"ok": False, "note": why}
            ctx = J["ctx"]
            k = ctx["seg_allele"]
            seg_in_cons = J["model"] == "SEG" and loc.get("consensus_allele", "") == L.alleles[k].name
            jr, n_el, elsewhere = fetch_any(bam, J["regions"], J["W"], args.max_reads)
            jal = [mappy.Aligner(seq=a.seq, preset="map-ont", best_n=5) for a in J["alleles"]]
            g = genotype_junctions(J, jr, jal, args.min_identity, threads=args.threads)
            est = junction_estimates(J["model"], g["type_counts"], seg_in_cons)
            rdir = os.path.join(sd["outdir"], "reads")
            os.makedirs(rdir, exist_ok=True)
            for a, nms in zip(J["alleles"], g["names"]):
                with open(os.path.join(rdir, f"{L.id}.junction.{a.name.replace(':', '_')}.txt"), "w") as fh:
                    fh.write("".join(x + "\n" for x in nms))
            size = ctx["e"] - ctx["s"]
            b = next((x for x in NULL_BINS if size <= x), NULL_BINS[-1])
            jkind = "INV" if J["model"] == "INV" else ("DEL" if seg_in_cons else "INS")
            floor = junction_floor((jkind, b))
            cnt = dict(zip(g["types"], g["type_counts"]))
            present, calls = [], []
            ma = args.min_alt_reads
            if J["model"] == "INV":
                v, lo, hi = est["inverted"]
                nref, ninv = cnt["REF_LEFT"] + cnt["REF_RIGHT"], cnt["INV_LEFT"] + cnt["INV_RIGHT"]
                if nref >= ma and 1 - hi > floor:
                    present.append(L.alleles[0].name)
                if ninv >= ma and lo > floor:
                    present.append(L.alleles[1].name)
                alt = est["inverted"]
                calls.append("inversion_mixture" if len(present) == 2 else "single_orientation" if present else "undetermined")
            else:
                xa, xlo, xhi = est["segment_absent"]
                npres, nabs, ne2s = (cnt["PRESENT_LEFT"] + cnt["PRESENT_RIGHT"]) / 2, cnt["ABSENT"], cnt["END_TO_START"]
                e2s = est.get("end_to_start_per_genome")
                absent_ok = nabs >= ma and xlo > floor
                present_ok = npres >= ma and 1 - xhi > floor
                e2s_ok = bool(e2s) and ne2s >= ma and e2s[1] > floor
                if "dup_allele" in ctx:
                    d = ctx["dup_allele"]
                    if e2s_ok:
                        present.append(L.alleles[d].name)
                    if e2s and present_ok and 1 - min(1.0, e2s[2]) > floor:
                        present.append(L.alleles[k].name)
                    alt = e2s if e2s else (None, 0.0, 1.0)
                else:
                    if present_ok:
                        present.append(L.alleles[k].name)
                    if absent_ok:
                        present.append(L.alleles[1 - k].name)
                    alt = est["segment_absent"] if k == 0 else est["segment_present"]
                if absent_ok and present_ok:
                    calls.append("segment_presence_mixture")
                elif absent_ok:
                    calls.append("segment_absent")
                elif present_ok:
                    calls.append("segment_present")
                if e2s_ok:
                    calls.append("end_to_start_junctions(tandem_copy_or_circular_form)")
                if absent_ok and e2s_ok:
                    calls.append("consistent_with_excision_and_circularisation")
            n = est["n"]
            if n < args.min_spanning:
                status = "insufficient_junction_reads"
            elif len(present) >= 2:
                status = "mixture"
            elif present:
                status = f"single:{present[0]}"
                if L.source == "discovered":
                    status = "consensus_discrepancy" if present[0] == L.alleles[1].name else "not_supported"
            else:
                status = "undetermined"
            fmt = lambda t: "NA" if t is None or t[0] is None else f"{t[0]:.4f} [{t[1]:.4f}-{min(t[2], 999):.4f}]"
            estxt = "; ".join(f"{kk}={fmt(vv)}" for kk, vv in est.items() if kk != "n")
            note = (f"junction window {J['W']} bp (anchors {J['RL']}/{J['RR']}); reads eligible {n_el}, used {len(jr)}, "
                    f"ambiguous {g['ambiguous']}, divergent {g['divergent']}, spanning none {g['unspanned']}"
                    + (f", primary_elsewhere {elsewhere}" if elsewhere else "")
                    + ("" if J["model"] == "INV" or seg_in_cons else "; end-to-start not assessable (segment not in consensus)"))
            allele_counts = dict(zip([a.name for a in J["alleles"]], g["counts"]))
            return {"ok": True, "status": status, "present": present, "counts": allele_counts, "estimates": est, "alt": alt,
                    "floor": floor, "call": ",".join(calls) or "none", "estimates_text": estxt, "note": note,
                    "fmt_alt": fmt(alt), "W": J["W"]}

        sres = []
        for L, loc in located:
            c = sd["cmap"][loc["contig"]]
            p1, p2 = (loc["start"] % c.len) + 1, ((loc["end"] - 1) % c.len) + 1
            ca = loc.get("consensus_allele", "")
            flags, present, parts, cov, genes = [], [], [], "", []
            n, n_el, reads, alle, win_note = 0, 0, [], {}, L.status_note
            floor = args.min_fraction
            g = None
            if L.W is None:
                status = "repeat_limited"
            else:
                regs = anchor_regions(L, loc)
                reads, n_el, fst = fetch_reads(bam, c, loc["start"] - 2000, loc["end"] + 2000, L.W, args.max_reads,
                                               regions=regs)
                if regs is None:
                    flags.append("consensus_matches_no_allele")
                elif fst["primary_elsewhere"]:
                    flags.append(f"primary_elsewhere:{fst['primary_elsewhere']}")
                al = [mappy.Aligner(seq=a.seq, preset="map-ont", best_n=5) for a in L.alleles]
                g = genotype(L, reads, al, args.min_identity, threads=args.threads)
                n = sum(g["counts"])
                floor = floors.get(sig(L, loc), args.min_fraction)
                if L.A_left > args.anchor + 100 or L.A_right > args.anchor + 100:
                    flags.append("anchor_extended_past_repeat")
                for k, a in enumerate(L.alleles):
                    lo, hi = wilson(g["counts"][k], n)
                    frac = g["counts"][k] / n if n else 0.0
                    parts.append(f"{a.name}={g['counts'][k]} ({100 * frac:.2f}% [{100 * lo:.2f}-{100 * hi:.2f}])")
                    if n and g["counts"][k] >= args.min_alt_reads and lo > floor:
                        present.append(a.name)
                    if g["counts"][k] >= 10 and min(g["fwd"][k], g["counts"][k] - g["fwd"][k]) < 0.1 * g["counts"][k]:
                        flags.append(f"strand_bias:{a.name}")
                    alle[a.name] = {"reads": g["counts"][k], "fraction": frac if n else None, "ci95": [lo, hi]}
                alt_ks = [k for k, a in enumerate(L.alleles) if a.name.startswith("ALT")]
                if L.source == "discovered" and len(alt_ks) > 1:
                    # alternative breakpoint hypotheses of one event: presence is tested on their sum
                    ta = sum(g["counts"][k] for k in alt_ks)
                    lo, hi = wilson(ta, n)
                    parts.append(f"ALT_total={ta} ({100 * ta / n if n else 0:.2f}% [{100 * lo:.2f}-{100 * hi:.2f}])")
                    alle["ALT_total"] = {"reads": ta, "fraction": ta / n if n else None, "ci95": [lo, hi]}
                    present = [x for x in present if not x.startswith("ALT")]
                    if n and ta >= args.min_alt_reads and lo > floor:
                        present.append("ALT")
                if n < args.min_spanning:
                    status = "insufficient_spanning_reads"
                elif len(present) >= 2:
                    status = "mixture"
                elif len(present) == 1:
                    status = f"single:{present[0]}"
                    if L.source == "discovered":
                        status = "consensus_discrepancy" if present[0] == "ALT" else "not_supported"
                else:
                    status = "undetermined"
                win_note = (f"window {L.W} bp (anchors {L.A_left}/{L.A_right}); reads eligible {n_el}, used {len(reads)}, "
                            f"ambiguous {g['ambiguous']}, divergent {g['divergent']}, spanning none {g['unspanned']}")
                if g["divergent"] and g["divergent"] >= 0.05 * max(1, n + g["divergent"]):
                    flags.append(f"divergent_reads:{g['divergent']}")
            if c.circular and (loc["start"] < 0 or loc["end"] > c.len):
                flags.append("spans_contig_origin")
            longest = max(L.alleles, key=lambda a: a.vlen)
            if L.kind in ("DEL", "INS", "COMPLEX", "DUP") and ca == longest.name and longest.vlen >= 1000:
                ss = loc["start"] + longest.vs if loc["strand"] == "+" else loc["end"] - longest.ve
                cr = coverage_ratio(bam, c, ss, ss + longest.vlen, sd["rmap"])
                cov = "" if cr is None else str(cr)
            segs = L.meta.get("segments") if L.source == "discovered" else None
            if segs and L.kind in ("DUP", "DEL") and ca == "REF":
                # explicit called segment in REF-allele coordinates (REF = consensus here)
                vs_, ve_ = segs[0]
                ss = loc["start"] + vs_ if loc["strand"] == "+" else loc["end"] - ve_
                if L.kind == "DUP":
                    cr = coverage_ratio(bam, c, ss, ss + (ve_ - vs_), sd["rmap"])
                    cov = "" if cr is None else str(cr)
                if feats:
                    genes = genes_in(feats, name_map.get(c.name, c.name), (ss % c.len) + 1, ((ss + ve_ - vs_ - 1) % c.len) + 1)
            elif feats and ca:
                aa = next(a for a in L.alleles if a.name == ca)
                if aa.vlen > 0:
                    ss = loc["start"] + aa.vs if loc["strand"] == "+" else loc["end"] - aa.ve
                    genes = genes_in(feats, name_map.get(c.name, c.name), (ss % c.len) + 1, ((ss + aa.vlen - 1) % c.len) + 1)
            # ---- where the event sits in this consensus, and the genes at its ends (annotation only)
            seg_txt, bp_txt, seg0 = "", "", None
            if segs and L.kind in ("DUP", "DEL") and ca == "REF":
                vs_, ve_ = segs[0]
                seg0 = (loc["start"] + vs_ if loc["strand"] == "+" else loc["end"] - ve_, ve_ - vs_)
            elif ca:
                aa_ = next((a for a in L.alleles if a.name == ca), None)
                if aa_ is not None:
                    seg0 = (loc["start"] + aa_.vs if loc["strand"] == "+" else loc["end"] - aa_.ve, aa_.vlen)
            if seg0 is not None:
                ss_, vl_ = seg0
                fid = name_map.get(c.name, c.name)
                if vl_ > 0:
                    S_, E_ = (ss_ % c.len) + 1, ((ss_ + vl_ - 1) % c.len) + 1
                    seg_txt = f"{loc['contig']}:{S_}-{E_}"
                    if feats:
                        bp_txt = segment_genes(feats, fid, loc["contig"], S_, E_, c.len, c.circular)
                else:
                    pj = ((ss_ - 1) % c.len) + 1 if c.circular else ss_
                    qj = (pj % c.len) + 1 if c.circular else ss_ + 1
                    seg_txt = f"{loc['contig']}:{pj}^{qj}"
                    if feats:
                        bp_txt = junction_genes(feats, fid, loc["contig"], pj, qj, c.len, c.circular)
            # ---- breakpoint-junction genotyping (events longer than the reads, or always on request)
            method = "spanning"
            jres = None
            # junction genotyping is primary when too few reads span the windows, or when the windows are so long
            # that the reads spanning them are too few to see minor versions (pre-specified rule, not data-driven)
            span_short = L.W is None or n < args.min_spanning or (L.W > args.max_span_window and args.junction_mode != "never")
            if args.junction_mode == "always" or (args.junction_mode == "auto" and span_short
                                                  and L.lmax >= 2 * args.anchor):
                jres = run_junction(L, loc)
                if jres["ok"] and span_short and not jres["status"].startswith("insufficient"):
                    method = "junction"
                    status, present, floor = jres["status"], jres["present"], jres["floor"]
            jcols = ["", "", "", "", "", ""]
            if jres is not None:
                if jres["ok"]:
                    jcols = ["; ".join(f"{kk}={vv}" for kk, vv in jres["counts"].items()), jres["fmt_alt"],
                             jres["estimates_text"], f"{jres['floor']:.4f}", jres["call"], jres["note"]]
                else:
                    jcols = ["", "", "", "", "not_applicable", jres["note"]]
            rows.append([s.name, L.id, L.source, L.kind, loc["contig"], p1, p2, L.lmax, L.dlen, loc["strand"], ca, n,
                         "; ".join(parts), f"{floor:.4f}", ",".join(present), cov, status, win_note, ",".join(flags),
                         "; ".join(genes), L.desc, method] + jcols + [seg_txt, bp_txt])
            rdir = os.path.join(sd["outdir"], "reads")
            os.makedirs(rdir, exist_ok=True)
            if g is not None:
                for k, a in enumerate(L.alleles):
                    with open(os.path.join(rdir, f"{L.id}.{a.name}.txt"), "w") as fh:
                        fh.write("".join(x + "\n" for x in g["names"][k]))
            jsum = None
            if jres is not None:
                jsum = {kk: vv for kk, vv in jres.items() if kk in ("ok", "status", "present", "counts", "floor", "call", "note")}
                if jres["ok"]:
                    jsum["estimates"] = {kk: (list(vv) if isinstance(vv, tuple) else vv) for kk, vv in jres["estimates"].items()}
            sres.append({"locus": L.id, "kind": L.kind, "source": L.source, "status": status, "method": method, "n": n,
                         "alleles": alle, "floor": floor, "coverage_ratio": cov, "genes": genes, "flags": flags,
                         "segment": seg_txt, "breakpoint_genes": bp_txt,
                         "contig": loc["contig"], "start": p1, "end": p2, "consensus_allele": ca,
                         "eligible_reads": n_el, "used_reads": len(reads), "window_bp": L.W, "junction": jsum})
        summary["samples"][s.name] = {"loci": sres, "discovery": sd["disc_note"], "merged_into": merged,
                                      "floors": {f"{k[0]}{k[1]}": v for k, v in floors.items()}}
        with open(os.path.join(sd["outdir"], "alleles.fasta"), "w") as fh:
            for L, loc in located:
                for a in L.alleles:
                    w = f" window={a.win[0] + 1}-{a.win[1]}" if a.win else ""
                    fh.write(f">{L.id}.{a.name} variable={a.vs + 1}-{a.ve}{w}\n{a.seq}\n")
        bam.close()

    hdr = ["sample", "locus", "source", "kind", "contig", "start", "end", "longest_variable_bp",
           "allele_length_difference_bp", "consensus_strand", "consensus_allele", "spanning_reads",
           "allele_reads_fraction_ci95", "detection_floor", "alleles_present", "coverage_ratio", "status", "window",
           "flags", "genes_in_consensus_segment", "description", "method", "junction_counts",
           "junction_alt_fraction_ci95", "junction_estimates", "junction_floor", "junction_call", "junction_note",
           "consensus_segment", "breakpoint_genes"]
    with open(os.path.join(args.outdir, "events.tsv"), "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(args.outdir, "null_loci.tsv"), "w") as fh:
        fh.write("\t".join(["sample", "signature", "null_locus", "position", "window_bp", "reads_eligible", "reads_used",
                            "consensus_reads", "fake_allele_reads", "ambiguous", "spanning_none", "fake_fraction"]) + "\n")
        for r in null_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    # one line per sample for pipeline logs / reports
    with open(os.path.join(args.outdir, "sample_summary.tsv"), "w") as fh:
        cats = ("mixture", "consensus_discrepancy", "single", "not_supported", "insufficient", "repeat_limited",
                "undetermined")
        fh.write("\t".join(["sample", "loci_genotyped", "known_loci", "discovered_loci"] + list(cats)
                           + ["reported_events"]) + "\n")
        for s in samples:
            rs = [r for r in rows if r[0] == s.name]
            cnt = {k: 0 for k in cats}
            rep = []
            for r in rs:
                st = r[16]
                key = ("single" if st.startswith("single:") else "insufficient" if st.startswith("insufficient")
                       else st if st in cnt else "undetermined")
                cnt[key] += 1
                if r[2] == "known" or st in ("mixture", "consensus_discrepancy"):
                    frac = r[23] if r[21] == "junction" else r[12]
                    rep.append(f"{r[1]} [{st}, {r[21]}] {frac}")
            fh.write("\t".join([s.name, str(len(rs)), str(sum(1 for r in rs if r[2] == "known")),
                                str(sum(1 for r in rs if r[2] == "discovered"))] + [str(cnt[k]) for k in cats]
                               + [" | ".join(rep) or "-"]) + "\n")
    summary["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(args.outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    log(f"done in {summary['elapsed_s']} s: {len(rows)} sample-locus rows")


if __name__ == "__main__":
    main()
