#!/usr/bin/env python3
"""polish_tools.py - helpers for consensus polishing (s12b_polish) and its validation.

Subcommands
  stats      one pass over a FASTQ: reads, bases and read groups (JSON)
  subsample  deterministic read-name-hash subsamples at target depths plus an optional
             disjoint hold-out set, all in one pass
  rgheader   @RG header lines (DS:basecall_model=...) for the read groups of a FASTQ
  changes    exact draft -> polished change table with read support from a BAM
  diff       small differences between two assemblies (minimap2 -cx asm5 --cs)
  inject     inject SNPs and homopolymer +-1 indels into a genome, with a truth table
  assess     score injected-error recovery from a diff of re-polished vs original
  rgmode     whether a FASTQ's reads carry RG tags: tags, none or mixed
  finalize   check a dorado polish result against its draft and change table, write the
             polished genome with the draft's FASTA headers (length= updated) and info.json
  install    replace consensus/<sample>.fasta with the polished genome (atomic) and record
             its checksum in info.json
  is-installed  exit 0 if the consensus file is the installed polished genome
  describe   one-line summary of a polished sample (for pipeline.log)

Read selection is by a keyed hash of the read name, so a read is either in or out of
a subsample regardless of file order, and subsamples at increasing depth are nested.
The hold-out set takes reads from the opposite end of the hash range, so it is
disjoint from every polishing subsample as long as max(depth) + holdout depth does
not exceed the total depth (checked).
"""
import argparse
import gzip
import hashlib
import json
import os
import random
import re
import subprocess
import sys

VERSION = '1.1'
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


# ----------------------------------------------------------------------------- io
def open_text(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode[0])


def read_fasta(path):
    recs, name, buf = [], None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if name is not None:
                    recs.append((name, ''.join(buf).upper()))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf).upper()))
    if not recs:
        raise SystemExit(f'no sequences in {path}')
    return recs


def write_fasta(path, recs, width=80):
    tmp = path + '.tmp'
    with open(tmp, 'w') as fh:
        for name, seq in recs:
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + '\n')
    os.replace(tmp, path)


def fastq_records(path):
    with open_text(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s, p, q = fh.readline(), fh.readline(), fh.readline()
            if not q or not h.startswith('@') or not p.startswith('+'):
                raise SystemExit(f'malformed or truncated FASTQ record in {path}: {h[:60]!r}')
            yield h, s, p, q


def hash_unit(seed, read_id):
    d = hashlib.blake2b(f'{seed}:{read_id}'.encode(), digest_size=8).digest()
    return int.from_bytes(d, 'big') / 18446744073709551616.0


def read_id(header):
    return header[1:].split(None, 1)[0]


RG_RE = re.compile(r'\tRG:Z:(\S+)')


# -------------------------------------------------------------------------- stats
def cmd_stats(a):
    n = bases = 0
    rgs = {}
    for h, s, _, _ in fastq_records(a.fastq):
        n += 1
        L = len(s) - 1
        bases += L
        m = RG_RE.search(h)
        key = m.group(1) if m else ''
        rgs[key] = rgs.get(key, 0) + 1
    out = {'fastq': os.path.abspath(a.fastq), 'reads': n, 'bases': bases, 'read_groups': rgs,
           'size_bytes': os.path.getsize(a.fastq), 'mtime': int(os.path.getmtime(a.fastq)), 'version': VERSION}
    with open(a.out, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: out[k] for k in ('reads', 'bases')}))


# ---------------------------------------------------------------------- subsample
def genome_size(path):
    return sum(len(s) for _, s in read_fasta(path))


def cmd_subsample(a):
    st = json.load(open(a.stats))
    if os.path.getsize(a.fastq) != st['size_bytes']:
        raise SystemExit('stats JSON does not describe this FASTQ (size differs); re-run stats')
    G = genome_size(a.draft)
    depths = [int(x) for x in a.depths.split(',') if x]
    total_depth = st['bases'] / G
    frac = {d: min(1.0, d * G / st['bases']) for d in depths}
    hq = min(1.0, a.holdout_depth * G / st['bases']) if a.holdout_depth else 0.0
    if hq and max(frac.values(), default=0) + hq > 1.0:
        raise SystemExit(f'hold-out ({a.holdout_depth}x) plus largest subsample ({max(depths)}x) exceeds '
                         f'the available depth ({total_depth:.0f}x); hold-out would overlap')
    outs = {d: open(f'{a.prefix}.d{d}.fastq', 'w') for d in depths}
    hold = open(f'{a.prefix}.holdout.fastq', 'w') if hq else None
    cnt = {d: [0, 0] for d in depths}
    hcnt = [0, 0]
    for h, s, p, q in fastq_records(a.fastq):
        u = hash_unit(a.seed, read_id(h))
        rec = h + s + p + q
        for d in depths:
            if u < frac[d]:
                outs[d].write(rec)
                cnt[d][0] += 1
                cnt[d][1] += len(s) - 1
        if hold is not None and u >= 1.0 - hq:
            hold.write(rec)
            hcnt[0] += 1
            hcnt[1] += len(s) - 1
    for fh in outs.values():
        fh.close()
    if hold is not None:
        hold.close()
    res = {'genome_size': G, 'total_depth': total_depth, 'seed': a.seed,
           'subsamples': {str(d): {'fraction': frac[d], 'reads': cnt[d][0], 'bases': cnt[d][1],
                                   'depth': cnt[d][1] / G, 'fastq': f'{a.prefix}.d{d}.fastq'} for d in depths}}
    if hold is not None:
        res['holdout'] = {'fraction': hq, 'reads': hcnt[0], 'bases': hcnt[1], 'depth': hcnt[1] / G,
                          'fastq': f'{a.prefix}.holdout.fastq'}
    with open(f'{a.prefix}.subsample.json', 'w') as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: round(v['depth'], 1) for k, v in res['subsamples'].items()} |
                     ({'holdout': round(res['holdout']['depth'], 1)} if hold is not None else {})))


# ----------------------------------------------------------------------- rgheader
MODEL_RE = re.compile(r'(dna_r[0-9.]+_e[0-9.]+_\d+bps_(?:fast|hac|sup)@v[0-9.]+)')
RUNID_RE = re.compile(r'^([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})_')


def cmd_rgheader(a):
    st = json.load(open(a.stats))
    lines = []
    for rg in sorted(st['read_groups']):
        if not rg:
            if not a.basecall_model:
                raise SystemExit('reads without an RG tag and no --basecall-model given')
            continue
        m = MODEL_RE.search(rg)
        model = m.group(1) if m else a.basecall_model
        if not model:
            raise SystemExit(f'cannot parse a basecall model from read group {rg!r}; pass --basecall-model')
        r = RUNID_RE.match(rg)
        ds = (f'runid={r.group(1)} ' if r else '') + f'basecall_model={model}'
        lines.append(f'@RG\tID:{rg}\tPL:ONT\tDS:{ds}\tSM:{a.sample}')
    if '' in st['read_groups']:
        lines.append(f'@RG\tID:{a.sample}_noRG\tPL:ONT\tDS:basecall_model={a.basecall_model}\tSM:{a.sample}')
    with open(a.out, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f'{len(lines)} read group(s)')


# ------------------------------------------------------------------------ changes
def cs_ops(cs):
    """Yield (op, value) from a short-form minimap2 cs string."""
    for m in re.finditer(r'(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})', cs):
        t = m.group(1)
        yield t[0], t[1:]


def contig_changes(draft, pol, threads=1):
    """Exact change list turning draft into pol, as (pos0, ref, alt) in draft coordinates."""
    if draft == pol:
        return []
    import mappy
    al = mappy.Aligner(seq=draft, preset='asm5', n_threads=threads)
    hits = [h for h in al.map(pol, cs=True) if h.is_primary and h.strand == 1]
    if len(hits) != 1:
        raise SystemExit(f'polished contig aligns to its draft in {len(hits)} forward primary pieces; '
                         'cannot derive an exact change list')
    h = hits[0]
    raw = []
    if h.r_st > 0 or h.q_st > 0:
        raw.append((0, draft[:h.r_st], pol[:h.q_st]))
    r, q = h.r_st, h.q_st
    for op, v in cs_ops(h.cs):
        if op == ':':
            n = int(v)
            r += n
            q += n
        elif op == '*':
            raw.append((r, v[0].upper(), v[1].upper()))
            r += 1
            q += 1
        elif op == '+':
            raw.append((r, '', v.upper()))
            q += len(v)
        elif op == '-':
            raw.append((r, v.upper(), ''))
            r += len(v)
        else:
            raise SystemExit('spliced cs operation in assembly alignment')
    if h.r_en < len(draft) or h.q_en < len(pol):
        raw.append((h.r_en, draft[h.r_en:], pol[h.q_en:]))
    # merge touching operations into one change
    merged = []
    for pos, ref, alt in sorted(raw, key=lambda x: (x[0], len(x[1]))):
        if merged and merged[-1][0] + len(merged[-1][1]) == pos:
            p0, r0, a0 = merged[-1]
            merged[-1] = (p0, r0 + ref, a0 + alt)
        else:
            merged.append((pos, ref, alt))
    # anchor pure insertions/deletions VCF-style
    out = []
    for pos, ref, alt in merged:
        if ref and alt:
            out.append((pos, ref, alt))
        elif pos > 0:
            out.append((pos - 1, draft[pos - 1] + ref, draft[pos - 1] + alt))
        else:
            out.append((0, ref + draft[len(ref)], alt + draft[len(ref)]))
    if apply_changes(draft, out) != pol:
        raise SystemExit('derived change list does not reproduce the polished sequence')
    return out


def apply_changes(seq, changes):
    s = seq
    for pos, ref, alt in sorted(changes, key=lambda x: x[0], reverse=True):
        if s[pos:pos + len(ref)] != ref:
            raise SystemExit(f'REF mismatch at {pos + 1}')
        s = s[:pos] + alt + s[pos + len(ref):]
    return s


def kind_of(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return 'SNP'
    if len(ref) == len(alt):
        return 'MNP'
    if len(alt) > len(ref) and alt.startswith(ref):
        return 'INS'
    if len(ref) > len(alt) and ref.startswith(alt):
        return 'DEL'
    return 'COMPLEX'


def hp_run(seq, i):
    """Length of the homopolymer run containing position i."""
    if i < 0 or i >= len(seq):
        return 0
    b = seq[i]
    s = e = i
    while s > 0 and seq[s - 1] == b:
        s -= 1
    while e + 1 < len(seq) and seq[e + 1] == b:
        e += 1
    return e - s + 1


def run_bounds(seq, i):
    b = seq[i]
    s = e = i
    while s > 0 and seq[s - 1] == b:
        s -= 1
    while e + 1 < len(seq) and seq[e + 1] == b:
        e += 1
    return s, e + 1


def query_window(read, ws, we):
    """Query slice of a read aligned over reference [ws, we), or None if not spanned."""
    if read.reference_start > ws or read.reference_end is None or read.reference_end < we:
        return None
    rpos, qpos = read.reference_start, 0
    q0 = q1 = None
    for op, ln in read.cigartuples:
        if op in (4,):                     # soft clip
            qpos += ln
            continue
        if op == 5:                        # hard clip
            continue
        if op in (0, 7, 8):                # M = X
            if q0 is None and rpos + ln > ws:
                q0 = qpos + max(0, ws - rpos)
            if q1 is None and rpos + ln >= we:
                q1 = qpos + (we - rpos)
                break
            rpos += ln
            qpos += ln
        elif op == 1:                      # I
            qpos += ln
        elif op in (2, 3):                 # D N
            if q0 is None and rpos + ln > ws:
                q0 = qpos
            if q1 is None and rpos + ln >= we:
                q1 = qpos
                break
            rpos += ln
    if q0 is None or q1 is None or q1 < q0:
        return None
    return read.query_sequence[q0:q1]


def vcf_changes(path, drafts, pols):
    """Change list from dorado's variants.vcf: every record (PASS and LowQual) is applied by
    dorado, so the list must reproduce the polished sequence exactly; None if it does not."""
    ch, meta = {}, {}
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        gt = f[9].split(':')[0] if len(f) > 9 else '1'
        idx = [int(x) for x in re.split(r'[/|]', gt) if x.isdigit()]
        if not idx or max(idx) == 0:
            continue
        alt = f[4].split(',')[max(idx) - 1].upper()
        ch.setdefault(f[0], []).append((int(f[1]) - 1, f[3].upper(), alt))
        meta[(f[0], int(f[1]) - 1)] = (f[5], f[6])
    for ctg in drafts:
        try:
            if apply_changes(drafts[ctg], ch.get(ctg, [])) != pols[ctg]:
                return None, None
        except SystemExit:
            return None, None
    return {c: sorted(v) for c, v in ch.items()}, meta


def cmd_changes(a):
    import pysam
    drafts = dict(read_fasta(a.draft))
    pols = dict(read_fasta(a.polished))
    if set(drafts) != set(pols):
        raise SystemExit('draft and polished contig names differ')
    bam = pysam.AlignmentFile(a.bam) if a.bam else None
    source = 'alignment'
    vch, meta = (vcf_changes(a.vcf, drafts, pols) if a.vcf and os.path.exists(a.vcf) else (None, None))
    if vch is not None:
        source = 'dorado_vcf'
    meta = meta or {}
    rows = []
    for ctg in drafts:
        d, p = drafts[ctg], pols[ctg]
        clist = vch.get(ctg, []) if vch is not None else contig_changes(d, p, a.threads)
        for pos, ref, alt in clist:
            # support window: change +- flank, widened to whole homopolymer runs at both ends
            ws, we = max(0, pos - a.flank), min(len(d), pos + len(ref) + a.flank)
            ws = run_bounds(d, ws)[0]
            we = run_bounds(d, we - 1)[1]
            wd = d[ws:we]
            wp = d[ws:pos] + alt + d[pos + len(ref):we]
            hp = max(hp_run(d, i) for i in range(pos, pos + max(1, len(ref))))
            if len(ref) != len(alt):
                hp = max(hp, hp_run(d, pos + 1) if pos + 1 < len(d) else 0)
            depth = nref = nalt = 0
            if bam is not None:
                for rd in bam.fetch(ctg, ws, we):
                    if rd.is_secondary or rd.is_unmapped or rd.mapping_quality < a.min_mapq:
                        continue
                    w = query_window(rd, ws, we)
                    if w is None:
                        continue
                    depth += 1
                    nref += w == wd
                    nalt += w == wp
            q, flt = meta.get((ctg, pos), ('.', '.'))
            rows.append((ctg, pos + 1, ref, alt, kind_of(ref, alt), hp, depth, nref, nalt, q, flt))
    with open(a.out, 'w') as fh:
        fh.write('contig\tpos\tref\talt\tkind\thp_len\tdepth\tref_reads\talt_reads\tqual\tfilter\n')
        for r in rows:
            fh.write('\t'.join(map(str, r)) + '\n')
    kinds = {}
    for r in rows:
        kinds[r[4]] = kinds.get(r[4], 0) + 1
    print(json.dumps({'changes': len(rows), 'kinds': kinds, 'source': source}))


# --------------------------------------------------------------------------- diff
DIFF_HEADER = 'q_contig\tq_pos\tr_contig\tr_pos\tkind\tlength\tq_allele\tr_allele\n'


def assembly_diffs(query, ref, threads=4, minimap2='minimap2'):
    """Small differences of query relative to ref from primary asm5 alignments."""
    res = subprocess.run([minimap2, '-c', '-x', 'asm5', '--cs', '--secondary=no', '-t', str(threads), ref, query],
                         capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit(f'minimap2 failed: {res.stderr[-500:]}')
    diffs, seen = [], set()
    aligned_q = 0
    for line in res.stdout.splitlines():
        f = line.split('\t')
        tags = dict(t.split(':', 2)[::2] for t in f[12:] if t.count(':') >= 2)
        if tags.get('tp') != 'P':
            continue
        qn, ql, qs, qe, strand, rn, rl, rs, re_ = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8])
        aligned_q += qe - qs
        cs = tags.get('cs', '')
        r = rs
        q = qs if strand == '+' else qe
        step = 1 if strand == '+' else -1
        for op, v in cs_ops(cs):
            if op == ':':
                n = int(v)
                r += n
                q += step * n
                continue
            if op == '*':
                qpos = q + 1 if strand == '+' else q
                row = (qn, qpos, rn, r + 1, 'SNP', 1, v[1].upper(), v[0].upper())
                r += 1
                q += step
            elif op == '+':
                qpos = q + 1 if strand == '+' else q
                row = (qn, qpos, rn, r, 'INS', len(v), v.upper(), '')
                q += step * len(v)
            elif op == '-':
                qpos = q + 1 if strand == '+' else q
                row = (qn, qpos, rn, r + 1, 'DEL', len(v), '', v.upper())
                r += len(v)
            else:
                continue
            key = (row[2], row[3], row[4], row[6], row[7])
            if key not in seen:
                seen.add(key)
                diffs.append(row)
    return diffs, aligned_q


def write_diffs(path, diffs):
    with open(path, 'w') as fh:
        fh.write(DIFF_HEADER)
        for r in diffs:
            fh.write('\t'.join(map(str, r)) + '\n')


def cmd_diff(a):
    diffs, aligned = assembly_diffs(a.query, a.ref, a.threads, a.minimap2)
    write_diffs(a.out, diffs)
    qlen = sum(len(s) for _, s in read_fasta(a.query))
    small = [d for d in diffs if d[5] <= 50]
    summ = {'small': len(small), 'snp': sum(d[4] == 'SNP' for d in small),
            'indel': sum(d[4] in ('INS', 'DEL') for d in small), 'large': len(diffs) - len(small),
            'query_aligned_fraction': aligned / qlen if qlen else 0.0}
    with open(a.out + '.json', 'w') as fh:
        json.dump(summ, fh)
    print(json.dumps(summ))


# ------------------------------------------------------------------------- inject
def cmd_inject(a):
    rng = random.Random(a.seed)
    recs = read_fasta(a.genome)
    total = sum(len(s) for _, s in recs)
    taken = []                                    # (contig, start, end) blocked intervals

    def free(ctg, s, e):
        return all(c != ctg or e + a.spacing <= s0 or s - a.spacing >= e0 for c, s0, e0 in taken)

    def pick_pos():
        x = rng.randrange(total)
        for name, seq in recs:
            if x < len(seq):
                return name, seq, x
            x -= len(seq)

    truth = []
    edits = {name: [] for name, _ in recs}
    tries = 0
    while sum(t[1] == 'SNP' for t in truth) < a.snps:
        tries += 1
        if tries > 10 ** 7:
            raise SystemExit('could not place SNPs')
        name, seq, i = pick_pos()
        if i < a.end_margin or i >= len(seq) - a.end_margin or hp_run(seq, i) > 1 or seq[i] not in 'ACGT':
            continue
        if hp_run(seq, i - 1) > 2 or hp_run(seq, i + 1) > 2 or not free(name, i, i + 1):
            continue
        alt = rng.choice([b for b in 'ACGT' if b != seq[i]])
        taken.append((name, i, i + 1))
        edits[name].append((i, seq[i], alt))
        truth.append([f'snp{len(truth) + 1}', 'SNP', name, i + 1, seq[i], alt])
    runs = []
    for name, seq in recs:
        i = a.end_margin
        while i < len(seq) - a.end_margin:
            s, e = run_bounds(seq, i)
            if e - s >= 3 and seq[i] in 'ACGT':
                runs.append((name, s, e))
            i = e
    rng.shuffle(runs)
    for name, s, e in runs:
        if sum(t[1] == 'HP' for t in truth) >= a.hps:
            break
        if not free(name, s, e):
            continue
        seq = dict(recs)[name]
        b = seq[s]
        if rng.random() < 0.5:
            edits[name].append((s, b, ''))            # delete one base of the run
            op = f'-1'
        else:
            edits[name].append((s, '', b))            # insert one base into the run
            op = f'+1'
        taken.append((name, s, e))
        truth.append([f'hp{sum(t[1] == "HP" for t in truth) + 1}', 'HP', name, s + 1, b * (e - s), op])
    if sum(t[1] == 'HP' for t in truth) < a.hps:
        raise SystemExit('not enough free homopolymer runs')
    out = []
    for name, seq in recs:
        s = seq
        for pos, ref, alt in sorted(edits[name], key=lambda x: x[0], reverse=True):
            s = s[:pos] + alt + s[pos + len(ref):]
        out.append((name, s))
    write_fasta(a.out_fasta, out)
    with open(a.out_truth, 'w') as fh:
        fh.write('id\tkind\tcontig\tpos\toriginal\tinjected\n')
        for t in truth:
            fh.write('\t'.join(map(str, t)) + '\n')
    print(json.dumps({'snps': a.snps, 'hps': a.hps, 'seed': a.seed}))


def cmd_assess(a):
    truth = [l.rstrip('\n').split('\t') for l in open(a.truth)][1:]
    diffs = [l.rstrip('\n').split('\t') for l in open(a.diff)][1:]   # query = re-polished, ref = original
    index = {}
    for t in truth:
        index.setdefault(t[2], []).append((int(t[3]), len(t[4]), t))
    unrestored, new = set(), []
    for d in diffs:
        rpos = int(d[3])
        hit = None
        for pos, ln, t in index.get(d[2], []):
            if pos - 10 <= rpos <= pos + ln + 10:
                hit = t
                break
        if hit is None:
            new.append(d)
        else:
            unrestored.add(hit[0])
    with open(a.out, 'w') as fh:
        fh.write('id\tkind\trestored\n')
        for t in truth:
            fh.write(f'{t[0]}\t{t[1]}\t{0 if t[0] in unrestored else 1}\n')
    with open(a.out_new, 'w') as fh:
        fh.write(DIFF_HEADER)
        for d in new:
            fh.write('\t'.join(d) + '\n')
    res = {k: sum(1 for t in truth if t[1] == k and t[0] not in unrestored) for k in ('SNP', 'HP')}
    res['new'] = len(new)
    print(json.dumps(res))



# ----------------------------------------------------------------- production steps
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_fasta_headers(path):
    """[(name, full header line without '>', seq)] in file order."""
    out, hdr, buf = [], None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if hdr is not None:
                    out.append((hdr.split()[0], hdr, ''.join(buf).upper()))
                hdr, buf = line[1:].strip(), []
            elif line.strip():
                buf.append(line.strip())
    if hdr is not None:
        out.append((hdr.split()[0], hdr, ''.join(buf).upper()))
    return out


def cmd_rgmode(a):
    st = json.load(open(a.stats))
    rgs = st['read_groups']
    tagged = sum(v for k, v in rgs.items() if k)
    untagged = rgs.get('', 0)
    print('tags' if tagged and not untagged else ('none' if untagged and not tagged else 'mixed'))


def cmd_finalize(a):
    draft = read_fasta_headers(a.draft)
    pol = dict(read_fasta(a.polished))
    names = [n for n, _, _ in draft]
    if len(set(names)) != len(names):
        raise SystemExit('duplicate contig names in the draft')
    if set(names) != set(pol):
        raise SystemExit(f'contig names differ: draft {sorted(names)[:5]} vs polished {sorted(pol)[:5]}')
    # the change table must turn the draft into the polished genome exactly
    ch = {}
    with open(a.changes) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {h: i for i, h in enumerate(hdr)}
        rows = [l.rstrip('\n').split('\t') for l in fh if l.strip()]
    for r in rows:
        ch.setdefault(r[ix['contig']], []).append((int(r[ix['pos']]) - 1, r[ix['ref']], r[ix['alt']]))
    for n, lst in ch.items():
        lst.sort()
        for (p0, r0, _), (p1, _, _) in zip(lst, lst[1:]):
            if p1 < p0 + len(r0):
                raise SystemExit(f'{n}: overlapping changes at {p0 + 1} and {p1 + 1}')
    for n, _, seq in draft:
        if apply_changes(seq, sorted(ch.get(n, []))) != pol[n]:
            raise SystemExit(f"{n}: applying the change table to the draft does not reproduce dorado's consensus")
    dlen = sum(len(sq) for _, _, sq in draft)
    plen = sum(len(sq) for sq in pol.values())
    rel = abs(plen - dlen) / max(1, dlen)
    if rel > a.max_length_change:
        raise SystemExit(f'polished genome length {plen} differs from the draft ({dlen}) by {rel:.5f} '
                         f'(> polish_max_length_change {a.max_length_change})')
    tmp = a.out_polished + '.tmp'
    with open(tmp, 'w') as fh:
        for n, h, _ in draft:
            seq = pol[n]
            h2 = re.sub(r'(^|\s)length=\d+', lambda m: f'{m.group(1)}length={len(seq)}', h)
            fh.write(f'>{h2}\n')
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + '\n')
    os.replace(tmp, a.out_polished)
    sub = json.load(open(a.subsample))
    st = json.load(open(a.stats))
    ss = sub['subsamples'][str(a.depth)]
    kinds, filt = {}, {}
    for r in rows:
        kinds[r[ix['kind']]] = kinds.get(r[ix['kind']], 0) + 1
        if 'filter' in ix:
            filt[r[ix['filter']]] = filt.get(r[ix['filter']], 0) + 1
    info = {'sample': a.sample, 'depth_target': a.depth, 'depth_realised': round(ss['depth'], 2),
            'reads_used': ss['reads'], 'bases_used': ss['bases'], 'reads_available': st['reads'],
            'seed': sub['seed'], 'draft_length': dlen, 'polished_length': plen,
            'changes': len(rows), 'change_kinds': kinds, 'change_filters': filt,
            'draft_sha256': sha256_file(a.draft), 'polished_sha256': sha256_file(a.out_polished),
            'dorado': a.dorado_version, 'polish_tools': VERSION}
    with open(a.out_info, 'w') as fh:
        json.dump(info, fh, indent=1)
    print(json.dumps({k: info[k] for k in ('changes', 'change_kinds', 'reads_used', 'depth_realised')}))


def cmd_install(a):
    info = json.load(open(a.info))
    if sha256_file(a.polished) != info['polished_sha256']:
        raise SystemExit('polished FASTA does not match info.json')
    tmp = a.consensus + '.polish_tmp'
    with open(a.polished, 'rb') as src, open(tmp, 'wb') as dst:
        dst.write(src.read())
    os.replace(tmp, a.consensus)
    info['installed_sha256'] = sha256_file(a.consensus)
    info['installed_as'] = os.path.abspath(a.consensus)
    with open(a.info + '.tmp', 'w') as fh:
        json.dump(info, fh, indent=1)
    os.replace(a.info + '.tmp', a.info)
    print('installed')


def cmd_describe(a):
    d = json.load(open(a.info))
    kinds = ', '.join(f'{k} {v}' for k, v in sorted(d['change_kinds'].items())) or 'none'
    filt = ', '.join(f'{k} {v}' for k, v in sorted(d.get('change_filters', {}).items()))
    print(f"{d['changes']} changes ({kinds}{'; dorado ' + filt if filt else ''}) from {d['reads_used']:,} reads "
          f"({d['depth_realised']}x); length {d['draft_length']:,} -> {d['polished_length']:,} bp")


def cmd_is_installed(a):
    """Exit 0 if the consensus is the installed polished genome or, with --as-draft, if the
    polishing results in --info were made from this consensus (review mode)."""
    try:
        info = json.load(open(a.info))
        key = 'draft_sha256' if a.as_draft else 'installed_sha256'
        ok = info.get(key) == sha256_file(a.consensus)
    except (OSError, ValueError, KeyError):
        ok = False
    sys.exit(0 if ok else 1)


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', action='version', version=VERSION)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('stats'); p.add_argument('--fastq', required=True); p.add_argument('--out', required=True)
    p = sp.add_parser('subsample')
    p.add_argument('--fastq', required=True); p.add_argument('--stats', required=True)
    p.add_argument('--draft', required=True); p.add_argument('--depths', required=True)
    p.add_argument('--holdout-depth', type=float, default=0); p.add_argument('--seed', default='polish1')
    p.add_argument('--prefix', required=True)
    p = sp.add_parser('rgheader')
    p.add_argument('--stats', required=True); p.add_argument('--sample', required=True)
    p.add_argument('--basecall-model', default=''); p.add_argument('--out', required=True)
    p = sp.add_parser('changes')
    p.add_argument('--draft', required=True); p.add_argument('--polished', required=True)
    p.add_argument('--bam', default=''); p.add_argument('--out', required=True)
    p.add_argument('--vcf', default='', help="dorado polish variants.vcf (preferred source of the change list)")
    p.add_argument('--flank', type=int, default=5); p.add_argument('--min-mapq', type=int, default=10)
    p.add_argument('--threads', type=int, default=4)
    p = sp.add_parser('diff')
    p.add_argument('--query', required=True); p.add_argument('--ref', required=True)
    p.add_argument('--out', required=True); p.add_argument('--threads', type=int, default=4)
    p.add_argument('--minimap2', default='minimap2')
    p = sp.add_parser('inject')
    p.add_argument('--genome', required=True); p.add_argument('--out-fasta', required=True)
    p.add_argument('--out-truth', required=True); p.add_argument('--snps', type=int, default=300)
    p.add_argument('--hps', type=int, default=300); p.add_argument('--seed', type=int, default=7)
    p.add_argument('--spacing', type=int, default=200); p.add_argument('--end-margin', type=int, default=1000)
    p = sp.add_parser('assess')
    p.add_argument('--truth', required=True); p.add_argument('--diff', required=True)
    p.add_argument('--out', required=True); p.add_argument('--out-new', required=True)
    p = sp.add_parser('rgmode'); p.add_argument('--stats', required=True)
    p = sp.add_parser('finalize')
    p.add_argument('--sample', required=True); p.add_argument('--draft', required=True)
    p.add_argument('--polished', required=True); p.add_argument('--changes', required=True)
    p.add_argument('--stats', required=True); p.add_argument('--subsample', required=True)
    p.add_argument('--depth', type=int, required=True)
    p.add_argument('--max-length-change', type=float, default=0.001)
    p.add_argument('--out-polished', required=True); p.add_argument('--out-info', required=True)
    p.add_argument('--dorado-version', default='')
    p = sp.add_parser('install')
    p.add_argument('--info', required=True); p.add_argument('--polished', required=True)
    p.add_argument('--consensus', required=True)
    p = sp.add_parser('is-installed')
    p.add_argument('--info', required=True); p.add_argument('--consensus', required=True)
    p.add_argument('--as-draft', action='store_true')
    p = sp.add_parser('describe'); p.add_argument('--info', required=True)
    a = ap.parse_args()
    {'stats': cmd_stats, 'subsample': cmd_subsample, 'rgheader': cmd_rgheader, 'changes': cmd_changes,
     'diff': cmd_diff, 'inject': cmd_inject, 'assess': cmd_assess, 'rgmode': cmd_rgmode,
     'finalize': cmd_finalize, 'install': cmd_install, 'is-installed': cmd_is_installed,
     'describe': cmd_describe}[a.cmd](a)


if __name__ == '__main__':
    main()
