#!/usr/bin/env python3
"""snv_tools.py - small variants between samples and against references (s18c_snv).

Subcommands
  repeats    self-alignment repeat mask of a genome (BED; same parameters as s18b_sv)
  mutate     mutated reference + truth VCF for caller validation (reads of T mapped to R;
             truth = the R->T differences in R coordinates)
  af-filter  Clair3 diploid output -> haploid calls: a het becomes hom-ALT (its highest-AF
             ALT) when AF >= threshold, otherwise it is dropped; hom-ALT calls are kept
  pass       keep PASS records, one ALT per record, haploid or diploid GT style
  combine    intersection or union of two normalised call sets (exact allele match)
  liftover   transfer CDS annotation from a sample genome to a closely related reference
  annotate   calls (normalised VCF) -> TSV with strand-resolved read support, repeat and
             homopolymer flags, and gene effect (translation table 11)
  vcfdist-counts  NONE-threshold SNP/INDEL counts from a vcfdist summary
  plan       choose the comparisons of a run: self, cross-sample (near-identical consensus
             genomes) and external reference (close enough to the sample), from k-mer
             containment and asm5 alignment identity / aligned fraction
  run        one comparison end to end: read subsample, minimap2 lr:hq mapping with read
             groups, the configured caller(s), normalisation, annotation
  collect    run-level comparisons.tsv and calls.tsv

Coordinates in TSV output are 1-based. Read support for each call compares every read's
sequence over a window (the call +-5 bp, widened to whole homopolymer runs) with the
reference window and the ALT window: exact matches are counted per strand.
"""
import argparse
import gzip
import json
import os
import random
import re
import sys

VERSION = '1.2'
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CODON = {}
_b = 'TCAG'
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b):
            CODON[x + y + z] = _aa[16 * i + 4 * j + k]
STARTS = {'ATG', 'GTG', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA'}     # table 11 initiation codons


def rc(s):
    return s.translate(COMP)[::-1]


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_fasta(path):
    recs, name, buf, hdr = [], None, [], ''
    with open_text(path) as fh:
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


def write_fasta(path, recs, width=80):
    with open(path + '.tmp', 'w') as fh:
        for name, seq, *_ in recs:
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + '\n')
    os.replace(path + '.tmp', path)


def is_circular(header, length):
    h = header.lower()
    if re.search(r'circular=(true|yes|y|1)\b', h):
        return True
    if re.search(r'circular=(false|no|n|0)\b', h):
        return False
    return length > 1_000_000          # complete chromosome-scale contigs from Autocycler are circular


# ------------------------------------------------------------------------ repeats
def merge_intervals(ivs):
    out = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def self_repeats(seq, chunk=1000, step=500, min_len=300, min_ident=0.85):
    """Intervals of `seq` that align elsewhere in `seq` (either strand); as in s18b_sv."""
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


def repeat_map(recs):
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
    if path and os.path.exists(path):
        for line in open(path):
            f = line.split()
            if len(f) >= 3 and not line.startswith(('#', 'track')):
                out.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return {k: merge_intervals(v) for k, v in out.items()}


def in_ivs(ivs, pos0):
    lo, hi = 0, len(ivs)
    while lo < hi:
        mid = (lo + hi) // 2
        if ivs[mid][1] <= pos0:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(ivs) and ivs[lo][0] <= pos0 < ivs[lo][1]


def cmd_repeats(a):
    recs = read_fasta(a.genome)
    rm = repeat_map(recs)
    with open(a.out, 'w') as fh:
        for name, _, _ in recs:
            for s, e in rm.get(name, []):
                fh.write(f'{name}\t{s}\t{e}\n')
    tot = sum(e - s for v in rm.values() for s, e in v)
    print(json.dumps({'repeat_bp': tot, 'intervals': sum(len(v) for v in rm.values())}))


# ------------------------------------------------------------------------- mutate
def hp_bounds(seq, i):
    b = seq[i]
    s = e = i
    while s > 0 and seq[s - 1] == b:
        s -= 1
    while e + 1 < len(seq) and seq[e + 1] == b:
        e += 1
    return s, e + 1


def cmd_mutate(a):
    """Edit T into R. Reads (from T) mapped to R should reveal the R->T differences."""
    rng = random.Random(a.seed)
    recs = read_fasta(a.genome)
    rmask = load_bed(a.repeats)
    seqs = {n: s for n, s, _ in recs}
    order = [n for n, _, _ in recs]
    total = sum(len(s) for s in seqs.values())
    blocked = {n: [] for n in order}

    def free(n, s, e):
        return all(e + a.spacing <= s0 or s - a.spacing >= e0 for s0, e0 in blocked[n])

    def rand_pos(in_repeat):
        for _ in range(200000):
            x = rng.randrange(total)
            for n in order:
                if x < len(seqs[n]):
                    break
                x -= len(seqs[n])
            if x < a.end_margin or x >= len(seqs[n]) - a.end_margin:
                continue
            if in_ivs(rmask.get(n, []), x) != in_repeat:
                continue
            return n, x
        raise SystemExit('could not place a variant (repeat fraction too small?)')

    edits = {n: [] for n in order}     # (t_start, t_end, r_seq) replace T[t_start:t_end] by r_seq
    n_rep_snp = round(a.snps * a.repeat_frac)
    n_rep_ind = round(a.indels * a.repeat_frac)
    plan = ([('SNP', True)] * n_rep_snp + [('SNP', False)] * (a.snps - n_rep_snp) +
            [('IND', True)] * n_rep_ind + [('IND', False)] * (a.indels - n_rep_ind))
    k_ind = 0
    for kind, rep in plan:
        for _ in range(100000):
            n, x = rand_pos(rep)
            t = seqs[n]
            if t[x] not in 'ACGT':
                continue
            if kind == 'SNP':
                if not free(n, x, x + 1):
                    continue
                edits[n].append((x, x + 1, rng.choice([b for b in 'ACGT' if b != t[x]]), 'SNP'))
                blocked[n].append((x, x + 1))
                break
            hp = (k_ind % 2 == 0)          # alternate: homopolymer-length change, random indel
            if hp:
                s, e = hp_bounds(t, x)
                if e - s < 3 or not free(n, s, e):
                    continue
                d = rng.choice([-2, -1, 1, 2])
                if e - s + d < 1:
                    continue
                rseq = t[s] * (e - s + d)
                edits[n].append((s, e, rseq, 'HP'))
                blocked[n].append((s, e))
            else:
                L = rng.randint(1, 20)
                if rng.random() < 0.5:           # R lacks T[x:x+L]  -> truth insertion
                    if not free(n, x, x + L):
                        continue
                    edits[n].append((x, x + L, '', 'IND'))
                    blocked[n].append((x, x + L))
                else:                            # R carries extra bases -> truth deletion
                    if not free(n, x, x + 1):
                        continue
                    ins = ''.join(rng.choice('ACGT') for _ in range(L))
                    edits[n].append((x, x, ins, 'IND'))
                    blocked[n].append((x, x + 1))
            k_ind += 1
            break
        else:
            raise SystemExit('variant placement failed')
    out_recs, truth = [], []
    for n in order:
        t = seqs[n]
        ed = sorted(edits[n])
        parts, last, rpos = [], 0, 0
        for ts, te, rseq, kind in ed:
            parts.append(t[last:ts])
            rpos += ts - last
            # truth record in R coordinates: R allele = rseq, T allele = t[ts:te]
            r_al, t_al = rseq, t[ts:te]
            if len(r_al) == len(t_al) == 1:
                truth.append((n, rpos + 1, r_al, t_al, kind))
            else:
                # anchor on the preceding base (identical in R and T)
                anc = t[ts - 1]
                truth.append((n, rpos, anc + r_al, anc + t_al, kind))
            parts.append(rseq)
            rpos += len(rseq)
            last = te
        parts.append(t[last:])
        out_recs.append((n, ''.join(parts)))
    write_fasta(a.out_ref, out_recs)
    with open(a.out_truth, 'w') as fh:
        fh.write('##fileformat=VCFv4.2\n')
        for n, s in out_recs:
            fh.write(f'##contig=<ID={n},length={len(s)}>\n')
        fh.write('##INFO=<ID=KIND,Number=1,Type=String,Description="SNP, HP (homopolymer length) or IND">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttruth\n')
        for n, p, r_al, t_al, kind in sorted(truth, key=lambda x: (order.index(x[0]), x[1])):
            fh.write(f'{n}\t{p}\t.\t{r_al}\t{t_al}\t60\tPASS\tKIND={kind}\tGT\t{a.gt}\n')
    with open(a.out_bed, 'w') as fh:
        for n, s in out_recs:
            fh.write(f'{n}\t0\t{len(s)}\n')
    rep_bed = load_bed(a.repeats)
    print(json.dumps({'truth': len(truth), 'snps': sum(t[4] == 'SNP' for t in truth),
                      'hp': sum(t[4] == 'HP' for t in truth), 'ind': sum(t[4] == 'IND' for t in truth),
                      'repeat_bp': sum(e - s for v in rep_bed.values() for s, e in v)}))


# --------------------------------------------------------------------- VCF tools
def vcf_iter(path):
    with open_text(path) as fh:
        for line in fh:
            yield line


def add_info(info, item):
    if not item:
        return info
    return item if info in ('.', '') else info + ';' + item


def parse_sample(fmt, val):
    return dict(zip(fmt.split(':'), val.split(':')))


def cmd_af_filter(a):
    out = open(a.out, 'w')
    kept = dropped = 0
    for line in vcf_iter(a.vcf):
        if line.startswith('#CHROM'):
            out.write('##INFO=<ID=C3AF,Number=1,Type=Float,Description="Clair3 allele fraction of the ALT">\n')
        if line.startswith('#'):
            out.write(line)
            continue
        f = line.rstrip('\n').split('\t')
        if f[6] not in ('PASS', '.'):
            continue
        smp = parse_sample(f[8], f[9])
        gt = re.split(r'[/|]', smp.get('GT', '.'))
        alts = f[4].split(',')
        try:
            afs = [float(x) for x in smp.get('AF', '').split(',')]
        except ValueError:
            afs = []
        if gt in (['1', '1'], ['1']):
            f[4], f[9], f[8] = alts[0], a.gt, 'GT'
            f[7] = add_info(f[7], f'C3AF={afs[0]:.4f}' if afs else '')
            out.write('\t'.join(f[:10]) + '\n')
            kept += 1
            continue
        gi = sorted({int(x) for x in gt if x.isdigit() and int(x) > 0})
        if not gi or not afs:
            dropped += 1
            continue
        # het (0/1, 1/2, 2/2 ...): choose the ALT with the highest AF among those in the genotype
        best = max(gi, key=lambda k: afs[k - 1] if k - 1 < len(afs) else -1)
        if best - 1 < len(afs) and afs[best - 1] >= a.af:
            f[7] = add_info(f[7], f'C3AF={afs[best - 1]:.4f}')
            f[4], f[9], f[8] = alts[best - 1], a.gt, 'GT'
            out.write('\t'.join(f[:10]) + '\n')
            kept += 1
        else:
            dropped += 1
    out.close()
    print(json.dumps({'kept': kept, 'dropped': dropped, 'af': a.af}))


def cmd_pass(a):
    out = open(a.out, 'w')
    n = 0
    seen_hdr = []
    for line in vcf_iter(a.vcf):
        if line.startswith('##FORMAT') and 'ID=GT' not in line:
            continue
        if line.startswith('#CHROM') and '##INFO=<ID=C3AF' not in ''.join(seen_hdr):
            out.write('##INFO=<ID=C3AF,Number=1,Type=Float,Description="Clair3 allele fraction of the ALT">\n')
        if line.startswith('#'):
            seen_hdr.append(line)
            out.write(line)
            continue
        f = line.rstrip('\n').split('\t')
        if f[6] not in ('PASS', '.'):
            continue
        if len(f) > 9:
            gt = parse_sample(f[8], f[9]).get('GT', '1')
            if gt in ('0', '0/0', '0|0', '.', './.'):
                continue
            idx = max(int(x) for x in re.split(r'[/|]', gt) if x.isdigit())
            afs = parse_sample(f[8], f[9]).get('AF', '').split(',')
            if afs and afs[0] not in ('', '.') and idx - 1 < len(afs):
                f[7] = add_info(f[7], f'C3AF={float(afs[idx - 1]):.4f}')
            f[4] = f[4].split(',')[idx - 1]
        f = f[:8] + ['GT', a.gt]
        out.write('\t'.join(f) + '\n')
        n += 1
    out.close()
    print(json.dumps({'records': n}))


def read_calls(path):
    hdr, recs = [], []
    for line in vcf_iter(path):
        if line.startswith('#'):
            hdr.append(line)
            continue
        f = line.rstrip('\n').split('\t')
        for alt in f[4].split(','):
            recs.append(f[:4] + [alt] + f[5:])
    return hdr, recs


def cmd_combine(a):
    ha, ra = read_calls(a.a)
    hb, rb = read_calls(a.b)
    ka = {(r[0], r[1], r[3], r[4]): r for r in ra}
    kb = {(r[0], r[1], r[3], r[4]): r for r in rb}
    if a.mode == 'int':
        keys = [k for k in ka if k in kb]
    else:
        keys = list(ka) + [k for k in kb if k not in ka]
    order = {}
    for line in ha:
        m = re.match(r'##contig=<ID=([^,>]+)', line)
        if m:
            order[m.group(1)] = len(order)
    keys.sort(key=lambda k: (order.get(k[0], 1 << 30), k[0], int(k[1])))
    with open(a.out, 'w') as fh:
        for line in ha:
            if line.startswith('#CHROM'):
                fh.write('##INFO=<ID=CALLERS,Number=1,Type=String,Description="Callers reporting the allele">\n')
                if not any(l.startswith('##INFO=<ID=C3AF') for l in ha):
                    fh.write('##INFO=<ID=C3AF,Number=1,Type=Float,Description="Clair3 allele fraction of the ALT">\n')
            fh.write(line)
        for k in keys:
            r = list(ka.get(k) or kb[k])
            info = r[7] if len(r) > 7 else '.'
            src = ('clair3,dorado' if k in ka and k in kb else ('clair3' if k in ka else 'dorado'))
            r[7] = (info + ';' if info not in ('.', '') else '') + f'CALLERS={src}'
            r = r[:8] + ['GT', a.gt]
            fh.write('\t'.join(r) + '\n')
    print(json.dumps({'mode': a.mode, 'a': len(ka), 'b': len(kb), 'out': len(keys)}))


def cmd_vcfdist_counts(a):
    rows = [l.rstrip('\n').split('\t') for l in open(a.summary) if l.strip()]
    hdr = rows[0]
    ix = {h: i for i, h in enumerate(hdr)}
    out = {}
    for r in rows[1:]:
        if r[ix['THRESHOLD']] == 'NONE' and r[ix['VAR_TYPE']] in ('SNP', 'INDEL'):
            out[r[ix['VAR_TYPE']]] = {k: int(r[ix[c]]) for k, c in (('truth_tp', 'TRUTH_TP'), ('query_tp', 'QUERY_TP'),
                                                                        ('fp', 'QUERY_FP'), ('fn', 'TRUTH_FN'))}
    print(json.dumps(out))


# ----------------------------------------------------------------- annotation
def gff_name_map(ref_recs, fna_path):
    """Annotation contig name -> FASTA contig name, by identical sequence (Bakta renames contigs)."""
    import hashlib
    if not fna_path or not os.path.exists(fna_path):
        return {}
    by_seq = {hashlib.sha1(seq.encode()).hexdigest(): name for name, seq, *_ in ref_recs}
    out = {}
    for n, sq, _ in read_fasta(fna_path):
        h = hashlib.sha1(sq.encode()).hexdigest()
        if h in by_seq:
            out[n] = by_seq[h]
    return out


def load_gff(path, name_map=None):
    """CDS/RNA features: contig -> sorted list of (start1, end1, strand, type, locus, gene, product)."""
    feats = {}
    name_map = name_map or {}
    if not path or not os.path.exists(path):
        return feats
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
            dec = lambda s: re.sub(r'%([0-9A-Fa-f]{2})', lambda m: chr(int(m.group(1), 16)), s)
            feats.setdefault(name_map.get(f[0], f[0]), []).append((int(f[3]), int(f[4]), f[6], f[2],
                                               dec(attr.get('locus_tag', attr.get('ID', ''))),
                                               dec(attr.get('gene', '')), dec(attr.get('product', ''))))
    for v in feats.values():
        v.sort()
    return feats


def effect(seq, feats, pos1, ref, alt):
    """Gene and effect for a VCF-style change at pos1 (1-based) on `seq`."""
    s0, e0 = pos1, pos1 + len(ref) - 1           # 1-based inclusive span of REF
    hits = [f for f in feats if f[0] <= e0 and f[1] >= s0]
    if not hits:
        near = min(feats, key=lambda f: min(abs(f[0] - s0), abs(f[1] - s0)), default=None)
        if near is None:
            return '', '', '', 'intergenic', ''
        d = min(abs(near[0] - s0), abs(near[1] - s0))
        return near[4], near[5], near[6], f'intergenic ({d} bp from {near[4] or near[5]})', ''
    f = hits[0]
    start, end, strand, ftype, locus, gene, product = f
    if ftype != 'CDS':
        return locus, gene, product, f'{ftype}', ''
    cds = seq[start - 1:end]
    if len(ref) != len(alt):
        diff = len(alt) - len(ref)
        kind = 'frameshift' if diff % 3 else ('inframe_insertion' if diff > 0 else 'inframe_deletion')
        return locus, gene, product, kind, ''
    new = seq[:pos1 - 1] + alt + seq[pos1 - 1 + len(ref):]
    ncds = new[start - 1:end]
    if strand == '-':
        cds, ncds = rc(cds), rc(ncds)
    aa0 = ''.join(CODON.get(cds[i:i + 3], 'X') for i in range(0, len(cds) - 2, 3))
    aa1 = ''.join(CODON.get(ncds[i:i + 3], 'X') for i in range(0, len(ncds) - 2, 3))
    if len(aa0) > 0 and cds[:3] in STARTS and aa0[0] != 'M':
        aa0 = 'M' + aa0[1:]
    if len(aa1) > 0 and ncds[:3] in STARTS and aa1[0] != 'M':
        aa1 = 'M' + aa1[1:]
    changes = [(i, x, y) for i, (x, y) in enumerate(zip(aa0, aa1)) if x != y]
    if not changes:
        if cds[:3] != ncds[:3]:
            return locus, gene, product, 'start_codon_changed', f'{cds[:3]}>{ncds[:3]}'
        return locus, gene, product, 'synonymous', ''
    i, x, y = changes[0]
    hgvs = ','.join(f'p.{p}{i + 1}{q}' for i, p, q in changes[:3])
    if i == 0 and cds[:3] in STARTS and ncds[:3] not in STARTS:
        return locus, gene, product, 'start_lost', hgvs
    if y == '*':
        return locus, gene, product, 'stop_gained', hgvs
    if x == '*':
        return locus, gene, product, 'stop_lost', hgvs
    return locus, gene, product, 'missense', hgvs


# ------------------------------------------------------------------- liftover
def cs_ops(cs):
    for m in re.finditer(r'(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})', cs):
        t = m.group(1)
        yield t[0], t[1:]


def cmd_liftover(a):
    """Sample CDS -> reference coordinates via asm5 alignment of the sample genome to the reference."""
    import mappy
    ref = read_fasta(a.ref)
    smp_recs = read_fasta(a.sample_fasta)
    smp = {n: s for n, s, _ in smp_recs}
    feats = load_gff(a.sample_gff, gff_name_map(smp_recs, a.sample_gff_fasta))
    out = []
    lifted = total = 0
    for rname, rseq, _ in ref:
        al = mappy.Aligner(seq=rseq, preset='asm5')
        for sname, sseq in smp.items():
            for h in al.map(sseq, cs=True):
                if not h.is_primary:
                    continue
                # build mapping sample pos0 -> ref pos0 along the alignment (reference orientation)
                mp = {}
                r, q = h.r_st, (h.q_st if h.strand == 1 else h.q_en - 1)
                step = 1 if h.strand == 1 else -1
                for op, v in cs_ops(h.cs):
                    if op == ':':
                        for _ in range(int(v)):
                            mp[q] = r
                            r += 1
                            q += step
                    elif op == '*':
                        mp[q] = r
                        r += 1
                        q += step
                    elif op == '+':
                        q += step * len(v)
                    elif op == '-':
                        r += len(v)
                for f in feats.get(sname, []):
                    st, en = f[0] - 1, f[1] - 1
                    if not (min(h.q_st, h.q_en) <= st and en < max(h.q_st, h.q_en)):
                        continue
                    total += 1
                    if st in mp and en in mp:
                        rs, re_ = sorted((mp[st], mp[en]))
                        strand = f[2] if h.strand == 1 else {'+': '-', '-': '+'}.get(f[2], f[2])
                        out.append((rname, rs + 1, re_ + 1, strand, f[3], f[4], f[5], f[6]))
                        lifted += 1
    with open(a.out, 'w') as fh:
        fh.write('##gff-version 3\n')
        for rname, s, e, strand, ftype, locus, gene, product in sorted(out):
            attr = f'locus_tag={locus}' + (f';gene={gene}' if gene else '') + (f';product={product}' if product else '')
            fh.write(f'{rname}\tliftover\t{ftype}\t{s}\t{e}\t.\t{strand}\t0\t{attr}\n')
    print(json.dumps({'features_in_aligned_blocks': total, 'lifted': lifted}))


# -------------------------------------------------------------------- annotate
def query_window(read, ws, we):
    if read.reference_start > ws or read.reference_end is None or read.reference_end < we:
        return None
    rpos, qpos = read.reference_start, 0
    q0 = q1 = None
    for op, ln in read.cigartuples:
        if op == 4:
            qpos += ln
            continue
        if op == 5:
            continue
        if op in (0, 7, 8):
            if q0 is None and rpos + ln > ws:
                q0 = qpos + max(0, ws - rpos)
            if q1 is None and rpos + ln >= we:
                q1 = qpos + (we - rpos)
                break
            rpos += ln
            qpos += ln
        elif op == 1:
            qpos += ln
        elif op in (2, 3):
            if q0 is None and rpos + ln > ws:
                q0 = qpos
            if q1 is None and rpos + ln >= we:
                q1 = qpos
                break
            rpos += ln
    if q0 is None or q1 is None or q1 < q0:
        return None
    return read.query_sequence[q0:q1]


def hp_len_at(seq, pos0, ref, alt):
    best = 0
    for i in range(max(0, pos0), min(len(seq), pos0 + max(1, len(ref)) + 1)):
        s, e = hp_bounds(seq, i)
        best = max(best, e - s)
    return best


def cmd_annotate(a):
    import pysam
    ref_recs = read_fasta(a.ref)
    seqs = {n: s for n, s, _ in ref_recs}
    feats = load_gff(a.gff, gff_name_map(ref_recs, a.gff_fasta))
    rmask = load_bed(a.repeats)
    bam = pysam.AlignmentFile(a.bam) if a.bam else None
    _, recs = read_calls(a.vcf)
    cols = ['comparison', 'contig', 'pos', 'ref', 'alt', 'vtype', 'qual', 'callers', 'clair3_af', 'depth',
            'ref_fwd', 'ref_rev', 'alt_fwd', 'alt_rev', 'in_repeat', 'hp_len', 'locus_tag', 'gene', 'product',
            'effect', 'aa_change', 'assembly_diff']
    import bisect
    adiffs, aligned, dpos = None, {}, {}
    if getattr(a, 'asm_query', ''):
        adiffs, aligned = assembly_diffs(a.asm_query, a.ref)
        for d in adiffs:
            dpos.setdefault(d[0], []).append(d[1])
    windows = []
    rows = []
    for r in recs:
        ctg, pos1, ref, alt = r[0], int(r[1]), r[3], r[4]
        seq = seqs[ctg]
        p0 = pos1 - 1
        if len(ref) == len(alt) == 1:
            vt = 'SNP'
        elif len(ref) == len(alt):
            vt = 'MNP'
        elif len(alt) > len(ref) and alt.startswith(ref[0]):
            vt = 'INS'
        elif len(ref) > len(alt) and ref.startswith(alt[0]):
            vt = 'DEL'
        else:
            vt = 'COMPLEX'
        info = dict(kv.split('=', 1) for kv in r[7].split(';') if '=' in kv) if len(r) > 7 else {}
        callers = info.get('CALLERS', a.callers)
        depth = rf = rr = af = ar = 0
        ws, we = max(0, p0 - 5), min(len(seq), p0 + len(ref) + 5)
        ws = hp_bounds(seq, ws)[0]
        we = hp_bounds(seq, we - 1)[1]
        windows.append((ctg, ws, we))
        if adiffs is None or not in_ivs(aligned.get(ctg, []), p0):
            ad = 'na'
        else:
            lst = dpos.get(ctg, [])
            k = bisect.bisect_left(lst, ws)
            ad = 'yes' if k < len(lst) and lst[k] <= we else 'no'
        if bam is not None:
            wref = seq[ws:we]
            walt = seq[ws:p0] + alt + seq[p0 + len(ref):we]
            for rd in bam.fetch(ctg, ws, we):
                if rd.is_secondary or rd.is_unmapped or rd.mapping_quality < a.min_mapq:
                    continue
                w = query_window(rd, ws, we)
                if w is None:
                    continue
                depth += 1
                if w == wref:
                    rr += rd.is_reverse
                    rf += not rd.is_reverse
                elif w == walt:
                    ar += rd.is_reverse
                    af += not rd.is_reverse
        if feats.get(ctg):
            locus, gene, product, eff, hgvs = effect(seq, feats[ctg], pos1, ref, alt)
        else:
            locus, gene, product, eff, hgvs = '', '', '', 'no_annotation', ''
        rows.append([a.comparison, ctg, pos1, ref, alt, vt, r[5], callers,
                     info.get('C3AF', ''), depth, rf, rr, af, ar,
                     'yes' if in_ivs(rmask.get(ctg, []), p0) else 'no', hp_len_at(seq, p0, ref, alt),
                     locus, gene, product, eff, hgvs, ad])
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for row in rows:
            fh.write('\t'.join(str(x).replace('\t', ' ') for x in row) + '\n')
    if adiffs is not None:
        wins = {}
        for ctg, ws, we in windows:
            wins.setdefault(ctg, []).append((ws, we))
        with open(os.path.join(os.path.dirname(os.path.abspath(a.out)), 'assembly_diffs.tsv'), 'w') as fh:
            fh.write('contig\tpos\tkind\tlength\tcalled\n')
            for ctg, pos0, kind, ln in adiffs:
                if ln > 50:
                    called = 'sv'
                else:
                    called = 'yes' if any(ws - 1 <= pos0 <= we for ws, we in wins.get(ctg, [])) else 'no'
                fh.write(f'{ctg}\t{pos0 + 1}\t{kind}\t{ln}\t{called}\n')
    print(json.dumps({'calls': len(rows), 'assembly_diffs': len(adiffs) if adiffs is not None else None}))


# ------------------------------------------------------------- comparison planning
def kmer_sketch(seq, k=21, scaled=200):
    """FracMinHash sketch: canonical k-mers kept when a fixed 64-bit mix of the k-mer falls
    below 2^64/scaled. Deterministic (no Python hash randomisation)."""
    import numpy as np
    lut = np.full(256, 255, dtype=np.uint8)
    for i, ch in enumerate(b'ACGT'):
        lut[ch] = i
    x = lut[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(x) - k + 1
    if n <= 0:
        return np.array([], dtype=np.uint64)
    bad = np.concatenate([[0], np.cumsum(x == 255)])
    okw = (bad[k:] - bad[:-k]) == 0
    xv = np.where(x == 255, 0, x).astype(np.uint64)
    fwd = np.zeros(n, dtype=np.uint64)
    rev = np.zeros(n, dtype=np.uint64)
    for j in range(k):
        fwd = (fwd << np.uint64(2)) | xv[j:j + n]
        rev = rev | ((np.uint64(3) - xv[j:j + n]) << np.uint64(2 * j))
    z = np.minimum(fwd, rev)[okw] + np.uint64(0x9E3779B97F4A7C15)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    z = z ^ (z >> np.uint64(31))
    return np.unique(z[z < np.uint64((1 << 64) // scaled)])


def containment(sa, sb):
    import numpy as np
    return float(len(np.intersect1d(sa, sb, assume_unique=True))) / max(1, len(sa))


def genome_similarity(query_fa, ref_fa, max_small=50):
    """asm5 primary alignments of query to reference. Identity counts only small differences
    (mismatches and indels <= max_small bp) over the aligned columns; larger indels are treated
    as unaligned sequence, so a strain with an extra prophage is not penalised in identity but
    in aligned fraction. Returns (identity, aligned fraction of query, aligned fraction of ref)."""
    import mappy
    al = mappy.Aligner(ref_fa, preset='asm5')
    if not al:
        raise SystemExit(f'cannot index {ref_fa}')
    q = read_fasta(query_fa)
    rlen = sum(len(s) for _, s, _ in read_fasta(ref_fa))
    qlen = sum(len(s) for _, s, _ in q)
    cols = diffs = big_q = big_r = 0
    qiv, riv = {}, {}
    for qname, qseq, _ in q:
        for h in al.map(qseq, cs=True):
            if not h.is_primary:
                continue
            for op, v in cs_ops(h.cs):
                if op == ':':
                    cols += int(v)
                elif op == '*':
                    cols += 1
                    diffs += 1
                elif len(v) <= max_small:
                    cols += len(v)
                    diffs += len(v)
                elif op == '+':
                    big_q += len(v)
                else:
                    big_r += len(v)
            qiv.setdefault(qname, []).append((h.q_st, h.q_en))
            riv.setdefault(h.ctg, []).append((h.r_st, h.r_en))
    aq = sum(e - s for v in qiv.values() for s, e in merge_intervals(v)) - big_q
    ar = sum(e - s for v in riv.values() for s, e in merge_intervals(v)) - big_r
    return (1 - diffs / cols if cols else 0.0), max(0, aq) / max(1, qlen), max(0, ar) / max(1, rlen)


def ref_name(path):
    b = os.path.basename(path)
    for ext in ('.gz', '.fasta', '.fna', '.fa', '.fas'):
        if b.endswith(ext):
            b = b[:-len(ext)]
    return b


def cmd_plan(a):
    """Comparisons to run: self (reads vs own consensus), cross (reads vs the consensus of a
    near-identical sample; groups by single linkage, the first sample in samplesheet order is
    the group's reference) and reference (reads vs an external genome close enough to the
    sample's consensus)."""
    samples = [l.rstrip('\n').split('\t') for l in open(a.samples) if l.strip()]
    refs = [l.rstrip('\n').split('\t') for l in open(a.references) if l.strip()] if a.references else []
    seqs = {s[0]: ''.join(x for _, x, _ in read_fasta(s[1])) for s in samples}
    sk = {k: kmer_sketch(v) for k, v in seqs.items()}
    rsk = {r[0]: kmer_sketch(''.join(x for _, x, _ in read_fasta(r[1]))) for r in refs}
    sim_rows, plan = [], []
    parent = {s[0]: s[0] for s in samples}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    by = {s[0]: s for s in samples}
    order = [s[0] for s in samples]
    for i, x in enumerate(order):
        for y in order[i + 1:]:
            cxy, cyx = containment(sk[x], sk[y]), containment(sk[y], sk[x])
            row = [x, y, 'sample', f'{cxy:.4f}', f'{cyx:.4f}', '', '', '', 'not_aligned (k-mer containment < 0.5)']
            if min(cxy, cyx) >= 0.5:
                ident, aq, ar = genome_similarity(by[y][1], by[x][1])
                ok = ident >= a.cross_min_identity and min(aq, ar) >= a.cross_min_aligned
                row[5:9] = [f'{ident:.6f}', f'{aq:.4f}', f'{ar:.4f}', 'cross' if ok else 'below_threshold']
                if ok and a.cross:
                    parent[find(y)] = find(x)
            sim_rows.append(row)
    for s in order:
        for r in refs:
            c = containment(sk[s], rsk[r[0]])
            row = [s, r[0], 'reference', f'{c:.4f}', '', '', '', '', 'not_aligned (k-mer containment < 0.3)']
            if c >= 0.3:
                ident, aq, ar = genome_similarity(by[s][1], r[1])
                ok = ident >= a.ref_min_identity and aq >= a.ref_min_aligned
                row[5:9] = [f'{ident:.6f}', f'{aq:.4f}', f'{ar:.4f}', 'reference' if ok else 'below_threshold']
                if ok:
                    gff = r[2] if len(r) > 2 else ''
                    smp = by[s]
                    plan.append([f'reference__{s}__{r[0]}', 'reference', s, r[0], r[1], smp[2], gff, '',
                                 '' if gff else smp[3], '' if gff else smp[1], '' if gff else smp[4], smp[1],
                                 f'{ident:.6f}', f'{aq:.4f}'])
            sim_rows.append(row)
    if a.self:
        for s in order:
            smp = by[s]
            plan.append([f'self__{s}', 'self', s, s, smp[1], smp[2], smp[3], smp[4], '', '', '', '', '1', '1'])
    groups = {}
    for s in order:
        groups.setdefault(find(s), []).append(s)
    for members in groups.values():
        rep = members[0]
        for m in members[1:]:
            sr = next(r for r in sim_rows if r[2] == 'sample' and {r[0], r[1]} == {rep, m})
            plan.append([f'cross__{m}__{rep}', 'cross', m, rep, by[rep][1], by[m][2], by[rep][3], by[rep][4],
                         '', '', '', by[m][1], sr[5] or '', sr[6] or ''])
    kind_order = {'self': 0, 'cross': 1, 'reference': 2}
    plan.sort(key=lambda p: (kind_order[p[1]], order.index(p[2]), p[3]))
    cols = ['comparison', 'kind', 'sample', 'reference', 'ref_fasta', 'reads', 'gff', 'gff_fasta', 'lift_gff',
            'lift_fasta', 'lift_gff_fasta', 'asm_query', 'identity', 'aligned_sample']
    with open(a.out, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for p in plan:
            fh.write('\t'.join(p) + '\n')
    with open(a.out_similarity, 'w') as fh:
        fh.write('query\ttarget\ttarget_kind\tkmer_containment_query\tkmer_containment_target\tidentity\t'
                 'aligned_fraction_query\taligned_fraction_target\tdecision\n')
        for r in sim_rows:
            fh.write('\t'.join(r) + '\n')
    print(json.dumps({k: sum(1 for p in plan if p[1] == k) for k in ('self', 'cross', 'reference')}))


# --------------------------------------------------------- assembly-level differences
def assembly_diffs(query_fa, ref_fa):
    """Differences of the query assembly relative to the reference from asm5 primary
    alignments (cs tag), in reference coordinates, plus the aligned reference intervals."""
    import mappy
    al = mappy.Aligner(ref_fa, preset='asm5')
    diffs, aligned = [], {}
    for _, qseq, _ in read_fasta(query_fa):
        for h in al.map(qseq, cs=True):
            if not h.is_primary:
                continue
            aligned.setdefault(h.ctg, []).append((h.r_st, h.r_en))
            r = h.r_st
            for op, v in cs_ops(h.cs):
                if op == ':':
                    r += int(v)
                elif op == '*':
                    diffs.append((h.ctg, r, 'SNP', 1))
                    r += 1
                elif op == '+':
                    diffs.append((h.ctg, r, 'INS', len(v)))
                elif op == '-':
                    diffs.append((h.ctg, r, 'DEL', len(v)))
                    r += len(v)
    return sorted(diffs), {k: merge_intervals(v) for k, v in aligned.items()}


# ---------------------------------------------------------------- run a comparison
def _sh(cmd, log, env=None):
    import subprocess
    with open(log, 'a') as lf:
        lf.write(f'$ {cmd}\n')
        lf.flush()
        r = subprocess.run(['bash', '-o', 'pipefail', '-c', cmd], stdout=lf, stderr=lf, env=env)
    if r.returncode:
        tail = open(log).read().splitlines()[-8:]
        raise SystemExit(f'command failed (rc={r.returncode}): {cmd[:200]}\n' + '\n'.join(tail))


def cmd_run(a):
    import shlex
    import shutil
    row = None
    with open(a.plan) as fh:
        cols = fh.readline().rstrip('\n').split('\t')
        for l in fh:
            f = dict(zip(cols, l.rstrip('\n').split('\t')))
            if f['comparison'] == a.comparison:
                row = f
    if row is None:
        raise SystemExit(f'{a.comparison} not in {a.plan}')
    od = os.path.abspath(a.outdir)
    os.makedirs(od, exist_ok=True)
    log = os.path.join(od, 'run.log')
    q = shlex.quote
    T = a.threads
    py = shlex.quote(sys.executable)
    pt = q(a.polish_tools)
    ref = os.path.join(od, 'ref.fasta')
    shutil.copyfile(row['ref_fasta'], ref)
    _sh(f'samtools faidx {q(ref)}', log)
    # reads: the sample's filtered reads, deterministic subsample at the target depth
    os.makedirs(a.stats_dir, exist_ok=True)
    stats = os.path.join(a.stats_dir, f"{row['sample']}.stats.json")
    fresh = False
    if os.path.exists(stats):
        st = json.load(open(stats))
        fresh = st.get('size_bytes') == os.path.getsize(row['reads']) and st.get('mtime') == int(os.path.getmtime(row['reads']))
    if not fresh:
        _sh(f"{py} {pt} stats --fastq {q(row['reads'])} --out {q(stats)}", log)
    _sh(f"{py} {pt} subsample --fastq {q(row['reads'])} --stats {q(stats)} --draft {q(ref)} --depths {a.depth} "
        f"--seed snv1 --prefix {q(od + '/sub')}", log)
    reads = f'{od}/sub.d{a.depth}.fastq'
    sub = json.load(open(f'{od}/sub.subsample.json'))
    import subprocess
    mode = subprocess.run([sys.executable, a.polish_tools, 'rgmode', '--stats', stats], capture_output=True,
                          text=True).stdout.strip()
    bam = f'{od}/aln.bam'
    if mode == 'tags':
        rg = f'{od}/rg.txt'
        _sh(f"{py} {pt} rgheader --stats {q(stats)} --sample {q(row['sample'])} --basecall-model {q(a.basecall_model)} "
            f"--out {q(rg)}", log)
        _sh(f'minimap2 -ax lr:hq -y -t {T} {q(ref)} {q(reads)} 2>> {q(od + "/minimap2.log")} '
            f'| samtools view -b -o {q(bam + ".u.bam")} -', log)
        _sh(f'{{ samtools view -H {q(bam + ".u.bam")} | awk \'!/^@RG/\'; cat {q(rg)}; }} > {q(bam + ".hdr")}', log)
        _sh(f'samtools reheader {q(bam + ".hdr")} {q(bam + ".u.bam")} | samtools sort -@ 4 -m 1G -o {q(bam)} -', log)
        os.remove(bam + '.u.bam')
        os.remove(bam + '.hdr')
    elif mode == 'none':
        if not a.basecall_model:
            raise SystemExit('reads carry no RG tags and no basecall model is configured')
        rgl = f"@RG\\tID:{row['sample']}\\tPL:ONT\\tDS:basecall_model={a.basecall_model}\\tSM:{row['sample']}"
        _sh(f"minimap2 -ax lr:hq -R {q(rgl)} -t {T} {q(ref)} {q(reads)} 2>> {q(od + '/minimap2.log')} "
            f"| samtools sort -@ 4 -m 1G -o {q(bam)} -", log)
    else:
        raise SystemExit('some reads carry RG tags and some do not; cannot build read groups for dorado')
    _sh(f'samtools index {q(bam)}', log)
    os.remove(reads)
    # callers (validated configurations; rule R2)
    cfgs = {'C3P': ('c3p',), 'C3AF': ('c3d',), 'DP': ('dp',), 'INT': ('c3d', 'dp'), 'UNI': ('c3d', 'dp')}
    if a.caller not in cfgs:
        raise SystemExit(f'unknown caller configuration {a.caller}')
    env = dict(os.environ)
    env['PATH'] = os.path.join(a.clair3_prefix, 'bin') + os.pathsep + env.get('PATH', '')
    env['CONDA_PREFIX'] = a.clair3_prefix
    norm = lambda i, o: _sh(f'bcftools norm -f {q(ref)} -m -any -c x -O v {q(i)} 2>/dev/null '
                             f'| bcftools sort -O v -o {q(o)} 2>/dev/null', log)
    made = {}
    for step in cfgs[a.caller]:
        if step in ('c3p', 'c3d'):
            cdir = f'{od}/clair3_{"haploid_precise" if step == "c3p" else "diploid"}'
            _sh(f'run_clair3.sh -b {q(bam)} -f {q(ref)} -m {q(a.clair3_model)} -t {T} -p ont -o {q(cdir)} '
                f'--include_all_ctgs --no_phasing_for_fa' + (' --haploid_precise' if step == 'c3p' else '')
                + f' > {q(cdir + ".log")} 2>&1', log, env=env)
            raw = f'{od}/{step}.vcf'
            _sh(f'zcat {q(cdir + "/merge_output.vcf.gz")} > {q(raw)}', log)
            shutil.rmtree(cdir, ignore_errors=True)
            filt = f'{od}/{step}.filtered.vcf'
            if step == 'c3p':
                cmd_pass(argparse.Namespace(vcf=raw, out=filt, gt='1'))
            else:
                cmd_af_filter(argparse.Namespace(vcf=raw, out=filt, af=a.af, gt='1'))
        else:
            ddir = f'{od}/dorado_polish'
            _sh(f'{q(a.dorado)} polish {q(bam)} {q(ref)} --bacteria --models-directory {q(a.models_dir)} -x {q(a.device)} '
                f'-t {T} -o {q(ddir)} --vcf --any-bam > {q(od + "/dorado_polish.log")} 2>&1', log)
            raw = f'{ddir}/variants.vcf'
            filt = f'{od}/dp.filtered.vcf'
            cmd_pass(argparse.Namespace(vcf=raw, out=filt, gt='1'))
        made[step] = f'{od}/{step}.norm.vcf'
        norm(filt, made[step])
    calls = f'{od}/calls.vcf'
    if a.caller in ('INT', 'UNI'):
        cmd_combine(argparse.Namespace(a=made['c3d'], b=made['dp'], mode=a.caller.lower()[:3], out=calls, gt='1'))
    else:
        shutil.copyfile(made[cfgs[a.caller][0]], calls)
    callers = 'dorado' if a.caller == 'DP' else 'clair3'
    # annotation: the reference's own GFF, or the sample's annotation lifted over to it
    rep = f'{od}/ref.repeats.bed'
    cmd_repeats(argparse.Namespace(genome=ref, out=rep))
    gff, gff_fa = row['gff'], row['gff_fasta']
    lift = {}
    if not gff and row['lift_gff']:
        gff, gff_fa = f'{od}/ref.liftover.gff3', ''
        cmd_liftover(argparse.Namespace(sample_fasta=row['lift_fasta'], sample_gff=row['lift_gff'],
                                        sample_gff_fasta=row['lift_gff_fasta'], ref=ref, out=gff))
        lift = {'liftover_from': row['lift_gff']}
    cmd_annotate(argparse.Namespace(vcf=calls, ref=ref, bam=bam, gff=gff, gff_fasta=gff_fa, repeats=rep,
                                    comparison=row['comparison'], callers=callers, min_mapq=10,
                                    asm_query=row['asm_query'], out=f'{od}/calls.tsv'))
    n = sum(1 for l in open(f'{od}/calls.tsv')) - 1
    info = {'comparison': row['comparison'], 'kind': row['kind'], 'sample': row['sample'],
            'reference': row['reference'], 'caller': a.caller, 'depth_target': a.depth,
            'depth_realised': round(sub['subsamples'][str(a.depth)]['depth'], 2),
            'reads_used': sub['subsamples'][str(a.depth)]['reads'], 'n_calls': n, 'annotation': gff or '',
            'snv_tools': VERSION} | lift
    with open(f'{od}/info.json', 'w') as fh:
        json.dump(info, fh, indent=1)
    open(f'{od}/.done', 'w').close()
    print(json.dumps({k: info[k] for k in ('comparison', 'n_calls', 'depth_realised')}))


def cmd_collect(a):
    """Run-level tables: comparisons.tsv (one row per planned comparison) and calls.tsv."""
    cols = None
    plan = []
    with open(a.plan) as fh:
        pc = fh.readline().rstrip('\n').split('\t')
        for l in fh:
            plan.append(dict(zip(pc, l.rstrip('\n').split('\t'))))
    out_rows, all_calls = [], []
    for p in plan:
        od = os.path.join(a.snv_dir, p['comparison'])
        ok = os.path.exists(os.path.join(od, '.done'))
        info = json.load(open(os.path.join(od, 'info.json'))) if ok else {}
        rows = []
        if ok:
            with open(os.path.join(od, 'calls.tsv')) as fh:
                hdr = fh.readline().rstrip('\n').split('\t')
                cols = cols or hdr
                rows = [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in fh if l.strip()]
        all_calls += rows
        coding = ('missense', 'stop_gained', 'stop_lost', 'start_lost', 'frameshift', 'inframe_insertion',
                  'inframe_deletion', 'start_codon_changed')
        adiff = {}
        if ok and os.path.exists(os.path.join(od, 'assembly_diffs.tsv')):
            with open(os.path.join(od, 'assembly_diffs.tsv')) as fh:
                h = fh.readline().rstrip('\n').split('\t')
                for l in fh:
                    d = dict(zip(h, l.rstrip('\n').split('\t')))
                    adiff[d['called']] = adiff.get(d['called'], 0) + 1
        na = 'NA'
        out_rows.append([p['comparison'], p['kind'], p['sample'], p['reference'],
                         info.get('depth_realised', na), len(rows),
                         sum(r['vtype'] == 'SNP' for r in rows), sum(r['vtype'] in ('INS', 'DEL') for r in rows),
                         sum(r['vtype'] in ('MNP', 'COMPLEX') for r in rows),
                         sum(r['effect'] in coding for r in rows),
                         sum(r.get('assembly_diff') == 'yes' for r in rows),
                         sum(r.get('assembly_diff') == 'no' for r in rows),
                         (adiff.get('yes', 0) + adiff.get('no', 0)) if p['kind'] != 'self' else na,
                         adiff.get('no', 0) if p['kind'] != 'self' else na,
                         p['identity'] or na, p['aligned_sample'] or na, info.get('caller', na),
                         'ok' if ok else 'failed'])
    with open(os.path.join(a.snv_dir, 'comparisons.tsv'), 'w') as fh:
        fh.write('comparison\tkind\tsample\treference\tdepth\tn_calls\tn_snp\tn_indel\tn_mnp_complex\tn_coding_change\t'
                 'calls_at_assembly_diff\tcalls_not_in_assembly_diff\tassembly_diffs\tassembly_diffs_without_call\t'
                 'identity\taligned_fraction_sample\tcaller\tstatus\n')
        for r in out_rows:
            fh.write('\t'.join(map(str, r)) + '\n')
    cols = cols or ['comparison', 'contig', 'pos', 'ref', 'alt', 'vtype', 'qual', 'callers', 'clair3_af', 'depth',
                    'ref_fwd', 'ref_rev', 'alt_fwd', 'alt_rev', 'in_repeat', 'hp_len', 'locus_tag', 'gene',
                    'product', 'effect', 'aa_change', 'assembly_diff']
    with open(os.path.join(a.snv_dir, 'calls.tsv'), 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in all_calls:
            fh.write('\t'.join(r.get(c, '') for c in cols) + '\n')
    print(json.dumps({'comparisons': len(out_rows), 'failed': sum(r[-1] == 'failed' for r in out_rows),
                      'calls': len(all_calls)}))


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', action='version', version=VERSION)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('repeats'); p.add_argument('--genome', required=True); p.add_argument('--out', required=True)
    p = sp.add_parser('mutate')
    p.add_argument('--genome', required=True); p.add_argument('--repeats', required=True)
    p.add_argument('--out-ref', required=True); p.add_argument('--out-truth', required=True)
    p.add_argument('--out-bed', required=True); p.add_argument('--snps', type=int, default=1500)
    p.add_argument('--indels', type=int, default=300); p.add_argument('--repeat-frac', type=float, default=0.10)
    p.add_argument('--spacing', type=int, default=50); p.add_argument('--end-margin', type=int, default=500)
    p.add_argument('--seed', type=int, default=11); p.add_argument('--gt', default='1')
    p = sp.add_parser('af-filter')
    p.add_argument('--vcf', required=True); p.add_argument('--out', required=True)
    p.add_argument('--af', type=float, default=0.65); p.add_argument('--gt', default='1')
    p = sp.add_parser('pass')
    p.add_argument('--vcf', required=True); p.add_argument('--out', required=True); p.add_argument('--gt', default='1')
    p = sp.add_parser('combine')
    p.add_argument('--a', required=True); p.add_argument('--b', required=True)
    p.add_argument('--mode', choices=['int', 'uni'], required=True); p.add_argument('--out', required=True)
    p.add_argument('--gt', default='1')
    p = sp.add_parser('vcfdist-counts'); p.add_argument('--summary', required=True)
    p = sp.add_parser('liftover')
    p.add_argument('--sample-fasta', required=True); p.add_argument('--sample-gff', required=True)
    p.add_argument('--ref', required=True); p.add_argument('--out', required=True)
    p.add_argument('--sample-gff-fasta', default='')
    p = sp.add_parser('annotate')
    p.add_argument('--vcf', required=True); p.add_argument('--ref', required=True)
    p.add_argument('--bam', default=''); p.add_argument('--gff', default='')
    p.add_argument('--repeats', default=''); p.add_argument('--comparison', required=True)
    p.add_argument('--gff-fasta', default='')
    p.add_argument('--callers', default='')
    p.add_argument('--asm-query', default='', help='sample assembly: flag calls at assembly-level differences')
    p.add_argument('--min-mapq', type=int, default=10); p.add_argument('--out', required=True)
    p = sp.add_parser('plan')
    p.add_argument('--samples', required=True, help='TSV: sample, consensus, reads, gff, annotation fna')
    p.add_argument('--references', default='', help='TSV: name, fasta, gff (optional)')
    p.add_argument('--self', action='store_true'); p.add_argument('--cross', action='store_true')
    p.add_argument('--cross-min-identity', type=float, default=0.999)
    p.add_argument('--cross-min-aligned', type=float, default=0.9)
    p.add_argument('--ref-min-identity', type=float, default=0.99)
    p.add_argument('--ref-min-aligned', type=float, default=0.8)
    p.add_argument('--out', required=True); p.add_argument('--out-similarity', required=True)
    p = sp.add_parser('run')
    p.add_argument('--plan', required=True); p.add_argument('--comparison', required=True)
    p.add_argument('--outdir', required=True); p.add_argument('--depth', type=int, required=True)
    p.add_argument('--caller', required=True); p.add_argument('--threads', type=int, default=8)
    p.add_argument('--af', type=float, default=0.65)
    p.add_argument('--polish-tools', required=True); p.add_argument('--stats-dir', required=True)
    p.add_argument('--basecall-model', default='')
    p.add_argument('--clair3-prefix', default=''); p.add_argument('--clair3-model', default='')
    p.add_argument('--dorado', default=''); p.add_argument('--models-dir', default='')
    p.add_argument('--device', default='cpu')
    p = sp.add_parser('collect')
    p.add_argument('--plan', required=True); p.add_argument('--snv-dir', required=True)
    a = ap.parse_args()
    {'repeats': cmd_repeats, 'mutate': cmd_mutate, 'af-filter': cmd_af_filter, 'pass': cmd_pass,
     'combine': cmd_combine, 'vcfdist-counts': cmd_vcfdist_counts, 'liftover': cmd_liftover,
     'annotate': cmd_annotate, 'plan': cmd_plan, 'run': cmd_run, 'collect': cmd_collect}[a.cmd](a)


if __name__ == '__main__':
    main()
