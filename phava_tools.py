#!/usr/bin/env python3
"""phava_tools.py - helpers for the s18e_phava stage (phase variation with PhaVa).

  gff        PhaVa-ready GFF3 from a Bakta or Prokka annotation: a region line for every
             consensus contig, every non-pseudo CDS with the consensus contig name and its
             coordinates unchanged, locus_tag and transl_table on every CDS and no FASTA
             section (PhaVa 0.2.3 reads only GFF files of this shape)
  subsample  read-name-hash subsample of a FASTQ to a maximum depth (the method of the
             polish and small-variant stages: a read is in or out by a keyed hash of its name)
  genotype   re-count every candidate with at least one PhaVa inverted read against its own
             two orientations only (PhaVa's references and counting rules, the same reads)
  annotate   per-sample inverton table: read counts with the reverse fraction, a Wilson 95 %
             interval, the call status and the gene context of every inverton
  collect    combined tables for all samples; invertons with the same segment sequence in
             several samples are linked

Why genotype: PhaVa maps all reads to the references of all candidates at once. Where
candidates are close their references share flank sequence, a read in the reference
orientation then fits two references equally (MAPQ 0) and is not counted, while a read with
the inverted segment fits only one; the inverted fraction is inflated (Cminuta01, 1:246611:
PhaVa 23/89, alone 24/809). Inverted reads are not lost this way, so candidates without a PhaVa
inverted read need no re-count; their counts stay PhaVa's (column counts_from).

PhaVa inverton identifiers are contig:a-b-c-d with 0-based starts and 1-based ends: left
inverted repeat a+1..b, invertible segment b+1..c, right inverted repeat c+1..d (1-based,
inclusive). PhaVa counts a read when its alignment to one orientation (and not the other)
runs from more than 20 bp before the left repeat to more than 20 bp after the right repeat.

Gene context. Features are the annotation's GFF lines except gene, region and source.
  class              the segment against the CDS (as PhaVa): intragenic (inside one CDS),
                     gene_boundary (overlaps a CDS partly) or intergenic
  genes_in_segment   every feature overlapping the segment
  breakpoint_genes   the features at the two ends of the segment, in the format of the SV
                     stage ([cut]: the feature extends beyond the segment)
  left/right gene    nearest gene (CDS, tRNA, rRNA, ncRNA, tmRNA) entirely outside the
                     inverton (repeats included) on each side, its strand and the number of
                     bases between it and the inverton; wraps on circular contigs
  promoter_targets   for intergenic invertons: the flanking genes whose 5' end faces the
                     inverton within --promoter-window bases (right gene on +, left gene on -),
                     i.e. the genes an invertible promoter in the segment could drive

PhaVa keeps candidates that overlap each other (one inversion can then be reported at more
than one candidate); overlapping_candidates lists them for every inverton.
"""
import argparse
import gzip
import hashlib
import json
import math
import os
import sys

VERSION = '1.1'
GENE_TYPES = ('CDS', 'tRNA', 'rRNA', 'ncRNA', 'tmRNA')
STATUSES = ('called', 'below_filter', 'no_reverse_reads', 'not_spanned')
CLASSES = ('intergenic', 'gene_boundary', 'intragenic')
RATIO_HEADER = ['Inverton', 'ForwardReadCount', 'ReverseReadCount', 'ProportionReverse', 'ReadFile',
                'Gene Overlaps', 'Upstream Gene', 'Upstream Strand', 'Distance to Upstream Gene',
                'Downstream Gene', 'Downstream Strand', 'Distance to Downstream Gene']
COLUMNS = ['sample', 'inverton', 'contig', 'start', 'end', 'left_ir', 'right_ir', 'segment', 'segment_bp',
           'ir_bp', 'forward_reads', 'reverse_reads', 'counted_reads', 'reverse_fraction', 'ci_low', 'ci_high',
           'status', 'counts_from', 'phava_forward_reads', 'phava_reverse_reads', 'class', 'genes_in_segment', 'breakpoint_genes', 'left_gene', 'left_strand',
           'left_distance', 'right_gene', 'right_strand', 'right_distance', 'promoter_targets',
           'phava_gene_overlaps', 'phava_upstream', 'phava_downstream', 'sequence_key', 'key_copies',
           'overlapping_candidates']
SUMMARY_COLUMNS = ['sample', 'candidates', 'spanned', 'with_reverse_reads', 'regenotyped', 'phava_called', 'called',
                   'called_intergenic',
                   'called_gene_boundary', 'called_intragenic', 'called_with_promoter_target',
                   'median_counted_reads', 'reads_used', 'depth_used', 'subsampled', 'min_reverse_reads',
                   'min_reverse_fraction', 'annotation', 'phava_version', 'tools_version']
COMBINED_EXTRA = ['in_other_samples']
GENOTYPE_COLUMNS = ['inverton', 'reads_tested', 'forward_reads', 'reverse_reads', 'phava_forward_reads',
                    'phava_reverse_reads']
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def rc(s):
    return s.translate(COMP)[::-1]


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_fasta(path):
    """[(name, sequence, header)] with the name the first word of the header."""
    recs, name, hdr, buf = [], None, '', []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    recs.append((name, ''.join(buf), hdr))
                hdr = line[1:].strip()
                name, buf = (hdr.split(None, 1)[0] if hdr else ''), []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf), hdr))
    return recs


def is_circular(header):
    toks = header.lower().split()
    return 'circular=true' in toks or 'topology=circular' in toks


def gff_fasta(path):
    """Sequences of the ##FASTA section of a GFF file."""
    recs, name, buf, on = [], None, [], False
    with open_text(path) as fh:
        for line in fh:
            if not on:
                on = line.startswith('##FASTA')
                continue
            if line.startswith('>'):
                if name is not None:
                    recs.append((name, ''.join(buf), ''))
                name, buf = line[1:].split(None, 1)[0], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, ''.join(buf), ''))
    return recs


def seq_md5(s):
    return hashlib.md5(s.upper().encode()).hexdigest()


def contig_map(cons, ann):
    """Consensus contig name -> annotation sequence id, matched by sequence. Without
    annotation sequences the names are taken to be the same."""
    if not ann:
        return {n: n for n, _, _ in cons}
    by = {}
    for n, s, _ in ann:
        by.setdefault(seq_md5(s), n)
    return {n: by[h] for n, s, _ in cons if (h := seq_md5(s)) in by}


def gff_rows(path):
    """(columns, attributes) of the GFF data lines before any ##FASTA section."""
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('##FASTA'):
                return
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            attrs = {}
            for kv in f[8].split(';'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    attrs.setdefault(k.strip(), v)
            yield f, attrs


def is_pseudo(attrs):
    return attrs.get('pseudo', '').lower() == 'true' or 'pseudogene' in attrs


def annotation_seqs(gff, fna):
    if fna and os.path.exists(fna):
        return read_fasta(fna)
    return gff_fasta(gff) if gff and os.path.exists(gff) else []


# ---------------------------------------------------------------------------- gff
def cmd_gff(a):
    cons = read_fasta(a.consensus)
    c2a = contig_map(cons, annotation_seqs(a.gff, a.fna))
    a2c = {v: k for k, v in c2a.items()}
    st = {'version': VERSION, 'gff': os.path.abspath(a.gff), 'consensus': os.path.abspath(a.consensus),
          'contigs': [{'name': n, 'length': len(s), 'annotation_seqid': c2a.get(n, '')} for n, s, _ in cons],
          'cds_in': 0, 'cds_out': 0, 'pseudo_skipped': 0, 'unmapped_skipped': 0, 'locus_tag_from_id': 0,
          'transl_table_added': 0}
    out = ['##gff-version 3']
    for n, s, _ in cons:
        out.append(f'{n}\tphava_tools\tregion\t1\t{len(s)}\t.\t+\t.\tID={n}')
    for f, attrs in gff_rows(a.gff):
        if f[2] != 'CDS':
            continue
        st['cds_in'] += 1
        if is_pseudo(attrs):
            st['pseudo_skipped'] += 1
            continue
        cn = a2c.get(f[0])
        if cn is None:
            st['unmapped_skipped'] += 1
            continue
        lt = attrs.get('locus_tag', '')
        if not lt:
            lt = attrs.get('ID', '') or f'cds{st["cds_out"] + 1}'
            st['locus_tag_from_id'] += 1
        tt = attrs.get('transl_table', '')
        if not tt:
            tt = str(a.transl_table)
            st['transl_table_added'] += 1
        # PhaVa takes the first attribute containing 'locus_tag' / 'transl_table', so the
        # attribute column carries only these two (the gene labels come from the annotation).
        out.append('\t'.join([cn] + f[1:8] + [f'ID={lt};locus_tag={lt};transl_table={tt}']))
        st['cds_out'] += 1
    tmp = a.out + '.tmp'
    with open(tmp, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    os.replace(tmp, a.out)
    with open(a.json or a.out + '.json', 'w') as fh:
        json.dump(st, fh, indent=1)
    print(json.dumps({k: st[k] for k in ('cds_in', 'cds_out', 'pseudo_skipped', 'unmapped_skipped',
                                          'locus_tag_from_id', 'transl_table_added')}))
    if st['cds_in'] and not st['cds_out']:
        print('warning: no CDS could be placed on a consensus contig; PhaVa will run without genes',
              file=sys.stderr)


# ---------------------------------------------------------------------- subsample
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


def cmd_subsample(a):
    G = sum(len(s) for _, s, _ in read_fasta(a.genome))
    n = bases = 0
    for _, s, _, _ in fastq_records(a.fastq):
        n += 1
        bases += len(s.rstrip('\n'))
    if not n:
        raise SystemExit(f'no reads in {a.fastq}')
    depth = bases / G
    res = {'version': VERSION, 'fastq': os.path.abspath(a.fastq), 'genome_size': G, 'reads': n, 'bases': bases,
           'depth': depth, 'max_depth': a.max_depth, 'seed': a.seed}
    if a.max_depth <= 0 or depth <= a.max_depth:
        res.update(subsampled=False, fraction=1.0, used_fastq=os.path.abspath(a.fastq), used_reads=n,
                   used_bases=bases, used_depth=depth)
    else:
        frac = a.max_depth * G / bases
        k = kb = 0
        tmp = a.out + '.tmp'
        with open(tmp, 'w') as out:
            for h, s, p, q in fastq_records(a.fastq):
                if hash_unit(a.seed, read_id(h)) < frac:
                    out.write(h + s + p + (q if q.endswith('\n') else q + '\n'))
                    k += 1
                    kb += len(s.rstrip('\n'))
        os.replace(tmp, a.out)
        res.update(subsampled=True, fraction=frac, used_fastq=os.path.abspath(a.out), used_reads=k,
                   used_bases=kb, used_depth=kb / G)
    with open(a.json, 'w') as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({'depth': round(depth, 1), 'used_depth': round(res['used_depth'], 1),
                      'subsampled': res['subsampled']}))


# ----------------------------------------------------------------------- genotype
def phava_count(sam, A, B, C, D, lo, max_mismatch):
    """PhaVa 0.2.3's counting (parseSam + computeRatios) for one inverton whose two references
    ('f', 'r') start at genome position lo: alignments with MAPQ > 2 and fewer mismatches plus
    deletions than max_mismatch of the segment and of the rest of the alignment; a read aligned
    to both orientations is dropped; the counted alignment spans the inverton by > 20 bp."""
    import pysam
    seg0 = B - lo                       # first segment base in the references
    seg1 = seg0 + (C - B)
    ist = A - lo                        # inverton start in the references
    ilen = D - A
    keep = {'f': {}, 'r': {}}
    with pysam.AlignmentFile(sam) as fh:
        for x in fh:
            if x.is_unmapped or x.mapping_quality <= 2:
                continue
            nm = nf = flen = 0
            for q, rp, base in x.get_aligned_pairs(with_seq=True):
                if rp is None:
                    continue
                if seg0 < rp < seg1:
                    if q is None or base.islower():
                        nm += 1
                else:
                    flen += 1
                    if q is None or base.islower():
                        nf += 1
            if nm < (seg1 - seg0) * max_mismatch and nf < flen * max_mismatch:
                keep[x.reference_name][x.query_name] = (x.reference_start, x.reference_end)
    n = {}
    for o, other in (('f', 'r'), ('r', 'f')):
        n[o] = sum(1 for name, (rs, re_) in keep[o].items()
                   if name not in keep[other] and rs < ist - 20 and re_ > ilen + ist + 20)
    return n['f'], n['r']


def cmd_genotype(a):
    import shutil
    import subprocess
    import tempfile
    import pysam
    mm2 = shutil.which('minimap2')
    if not mm2:
        raise SystemExit('minimap2 not on PATH')
    cons = {n: s for n, s, _ in read_fasta(a.consensus)}
    targets = [r for r in read_ratio(a.ratio) if int(r['ReverseReadCount']) >= 1]
    sub = json.load(open(a.subsample_json)) if a.subsample_json and os.path.exists(a.subsample_json) else {}
    frac, seed = float(sub.get('fraction', 1.0)), sub.get('seed', 'phava1')
    names = {}
    with pysam.AlignmentFile(a.bam) as bam:
        for r in targets:
            contig, A, B, C, D = parse_inverton(r['Inverton'])
            sel = set()
            # every alignment within the flanks: a read carrying the inversion can be split or
            # clipped at the inverton in the consensus mapping
            for x in bam.fetch(contig, max(0, A - a.flank), min(len(cons[contig]), D + a.flank)):
                if x.is_secondary or x.is_unmapped:
                    continue
                if frac < 1 and hash_unit(seed, x.query_name) >= frac:
                    continue
                sel.add(x.query_name)
            names[r['Inverton']] = sel
    wanted = set().union(*names.values()) if names else set()
    seqs = {}
    if wanted:
        for h, s, _, _ in fastq_records(a.fastq):
            rid = read_id(h)
            if rid in wanted:
                seqs[rid] = s.rstrip('\n')
    rows = []
    with tempfile.TemporaryDirectory(dir=a.tmp or None) as tmp:
        for r in targets:
            contig, A, B, C, D = parse_inverton(r['Inverton'])
            g = cons[contig]
            lo = max(0, A - a.flank)
            left, right = g[lo:A], g[D:D + a.flank]
            ref = os.path.join(tmp, 'ref.fa')
            with open(ref, 'w') as fh:
                fh.write(f'>f\n{left}{g[A:B]}{g[B:C]}{g[C:D]}{right}\n>r\n{left}{g[A:B]}{rc(g[B:C])}{g[C:D]}{right}\n')
            rd = os.path.join(tmp, 'reads.fa')
            use = sorted(n for n in names[r['Inverton']] if n in seqs)
            with open(rd, 'w') as fh:
                for n in use:
                    fh.write(f'>{n}\n{seqs[n]}\n')
            sam = os.path.join(tmp, 'aln.sam')
            if use:
                subprocess.run([mm2, '-a', '--MD', '-t', str(a.threads), '-o', sam, ref, rd], check=True,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                nf, nr = phava_count(sam, A, B, C, D, lo, a.max_mismatch)
            else:
                nf = nr = 0
            rows.append({'inverton': r['Inverton'], 'reads_tested': len(use), 'forward_reads': nf, 'reverse_reads': nr,
                         'phava_forward_reads': r['ForwardReadCount'], 'phava_reverse_reads': r['ReverseReadCount']})
    write_tsv(a.out, rows, GENOTYPE_COLUMNS)
    print(json.dumps({'regenotyped': len(rows), 'reads': len(seqs)}))

# ----------------------------------------------------------------------- annotate
def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, (c - h) / d), min(1.0, (c + h) / d)


def parse_inverton(inv):
    contig, coords = inv.rsplit(':', 1)
    a, b, c, d = (int(x) for x in coords.split('-'))
    if not (a < b <= c < d):
        raise SystemExit(f'unexpected inverton coordinates: {inv}')
    return contig, a, b, c, d


def read_ratio(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        if hdr != RATIO_HEADER:
            raise SystemExit(f'{path}: not a PhaVa 0.2.3 ratio table (header {hdr})')
        rows = []
        for line in fh:
            if line.strip():
                f = line.rstrip('\n').split('\t')
                f += [''] * (len(hdr) - len(f))
                rows.append(dict(zip(hdr, f)))
    return rows


def load_features(path):
    """seqid -> [(start, end, type, locus_tag, gene, product, strand, pseudo)], 1-based inclusive."""
    feats = {}
    if not path or not os.path.exists(path):
        return feats
    for f, attrs in gff_rows(path):
        if f[2] in ('gene', 'region', 'source'):
            continue
        prod = attrs.get('product', '').replace('%2C', ',').replace('%3B', ';')
        feats.setdefault(f[0], []).append((int(f[3]), int(f[4]), f[2], attrs.get('locus_tag', attrs.get('ID', '')),
                                           attrs.get('gene', ''), prod, f[6], is_pseudo(attrs)))
    return feats


def feat_label(f):
    """locus_tag(gene):product, the form used in every gene column of the pipeline."""
    return f"{f[3]}{'(' + f[4] + ')' if f[4] else ''}:{f[5]}"


def _covers(f, pos, clen):
    """Feature f contains 1-based pos; a feature with end > contig length crosses the origin."""
    return f[0] <= pos <= f[1] or (f[1] > clen and pos <= f[1] - clen)


def _cut(f, pos, clen, side):
    """Feature f (which contains pos) extends beyond a segment that starts or ends at pos."""
    if f[0] <= pos <= f[1]:
        return f[0] < pos if side == 'start' else f[1] > pos
    return True if side == 'start' else (f[1] - clen) > pos


def _overlaps(f, s, e, clen):
    """Feature f overlaps s..e (1-based, s <= e)."""
    if f[1] >= s and f[0] <= e:
        return True
    return f[1] > clen and s <= f[1] - clen


def _nearest(fl, lo, hi, clen, circular):
    """Nearest feature ending at or before lo and nearest starting at or after hi, with the
    number of bases between it and lo..hi (exclusive); wraps on circular contigs."""
    left = [f for f in fl if f[1] <= lo]
    right = [f for f in fl if f[0] >= hi]
    fL = fR = None
    gL = gR = None
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
    return fL, gL, fR, gR


def _flank_text(fl, lo, hi, clen, circular):
    fL, gL, fR, gR = _nearest(fl, lo, hi, clen, circular)
    lt = f'left {feat_label(fL)}, {gL} bp' if fL is not None else 'left none'
    rt = f'right {feat_label(fR)}, {gR} bp' if fR is not None else 'right none'
    return f'{lt}; {rt}'


def segment_genes(fl, ctg, S, E, clen, circular):
    """Features at both ends of segment S..E (1-based, inclusive), as in the SV stage."""
    def at(pos, side):
        inside = sorted((f for f in fl if _covers(f, pos, clen)), key=lambda f: (f[0], f[1], f[3]))
        if inside:
            return 'in ' + ', '.join(feat_label(f) + (' [cut]' if _cut(f, pos, clen, side) else '') for f in inside)
        return 'intergenic (' + _flank_text(fl, pos - 1, pos + 1, clen, circular) + ')'
    return f'start {ctg}:{S} {at(S, "start")} | end {ctg}:{E} {at(E, "end")}'


def classify(fl, S, E, clen):
    """Segment S..E against the non-pseudo CDS: intragenic, gene_boundary or intergenic."""
    cds = [f for f in fl if f[2] == 'CDS' and not f[7]]
    for f in cds:
        if f[0] <= S and E <= f[1]:
            return 'intragenic'
        if f[1] > clen and S <= E <= f[1] - clen:
            return 'intragenic'
    if any(_overlaps(f, S, E, clen) for f in cds):
        return 'gene_boundary'
    return 'intergenic'


def compact_overlaps(txt):
    out = []
    for item in txt.split(';'):
        if not item.strip():
            continue
        gene, _, kind = item.partition(',')
        if kind.startswith('intragenic,'):
            out.append(f'{gene} intragenic ({kind.split(",")[1]})')
        elif kind.startswith('partial overlap'):
            out.append(f'{gene} partial overlap ({kind.split(",")[-1].strip()})')
        else:
            out.append(kind if not gene else f'{gene} {kind}')
    return '; '.join(out)


def phava_neighbour(gene, strand, dist):
    if not gene or gene == 'NA' or dist in ('', '1000000000'):
        return ''
    return f'{gene} ({strand}, {dist} bp)'


def fmt(x, nd=4):
    return 'NA' if x is None else f'{x:.{nd}f}'


def annotate_rows(sample, ratio_rows, cons, feats, c2a, min_rev, min_frac, window, geno=None):
    seqs = {n: (s, is_circular(h)) for n, s, h in cons}
    order = {n: i for i, (n, _, _) in enumerate(cons)}
    rows = []
    for r in ratio_rows:
        contig, a, b, c, d = parse_inverton(r['Inverton'])
        if contig not in seqs:
            raise SystemExit(f'{r["Inverton"]}: contig not in the consensus')
        seq, circ = seqs[contig]
        clen = len(seq)
        pF, pR = int(r['ForwardReadCount']), int(r['ReverseReadCount'])
        g = (geno or {}).get(r['Inverton'])
        F, R = (int(g['forward_reads']), int(g['reverse_reads'])) if g else (pF, pR)
        n = F + R
        frac = R / n if n else None
        lo, hi = wilson(R, n) if n else (None, None)
        if n == 0:
            status = 'not_spanned'
        elif R == 0:
            status = 'no_reverse_reads'
        elif R >= min_rev and R / n >= min_frac:
            status = 'called'
        else:
            status = 'below_filter'
        S, E = b + 1, c
        m = seq[b:c].upper()
        row = {'sample': sample, 'inverton': r['Inverton'], 'contig': contig, 'start': a + 1, 'end': d,
               'left_ir': f'{a + 1}-{b}', 'right_ir': f'{c + 1}-{d}', 'segment': f'{S}-{E}', 'segment_bp': E - S + 1,
               'ir_bp': b - a, 'forward_reads': F, 'reverse_reads': R, 'counted_reads': n,
               'reverse_fraction': fmt(frac), 'ci_low': fmt(lo), 'ci_high': fmt(hi), 'status': status,
               'counts_from': 'regenotyped' if g else 'phava', 'phava_forward_reads': pF, 'phava_reverse_reads': pR,
               'phava_gene_overlaps': compact_overlaps(r['Gene Overlaps']),
               'phava_upstream': phava_neighbour(r['Upstream Gene'], r['Upstream Strand'], r['Distance to Upstream Gene']),
               'phava_downstream': phava_neighbour(r['Downstream Gene'], r['Downstream Strand'],
                                                   r['Distance to Downstream Gene']),
               'sequence_key': hashlib.sha1(min(m, rc(m)).encode()).hexdigest()[:16], '_order': (order[contig], a)}
        sid = c2a.get(contig)
        if sid is None or sid not in feats:
            row.update({k: 'NA' for k in ('class', 'genes_in_segment', 'breakpoint_genes', 'left_gene', 'left_strand',
                                          'left_distance', 'right_gene', 'right_strand', 'right_distance',
                                          'promoter_targets')})
        else:
            fl = feats[sid]
            genes = [f for f in fl if f[2] in GENE_TYPES]
            cls = classify(fl, S, E, clen)
            seg = sorted((f for f in fl if _overlaps(f, S, E, clen)), key=lambda f: (f[0], f[1], f[3]))
            fL, gL, fR, gR = _nearest(genes, a, d + 1, clen, circ)
            targets = []
            if cls == 'intergenic':
                if fL is not None and fL[6] == '-' and gL <= window:
                    targets.append(f'{feat_label(fL)} (-, {gL} bp)')
                if fR is not None and fR[6] == '+' and gR <= window:
                    targets.append(f'{feat_label(fR)} (+, {gR} bp)')
            row.update({'class': cls, 'genes_in_segment': '; '.join(feat_label(f) for f in seg),
                        'breakpoint_genes': segment_genes(fl, contig, S, E, clen, circ),
                        'left_gene': feat_label(fL) if fL else '', 'left_strand': fL[6] if fL else '',
                        'left_distance': gL if fL else '', 'right_gene': feat_label(fR) if fR else '',
                        'right_strand': fR[6] if fR else '', 'right_distance': gR if fR else '',
                        'promoter_targets': '; '.join(targets)})
        rows.append(row)
    rows.sort(key=lambda x: x['_order'])
    copies = {}
    for x in rows:
        copies[x['sequence_key']] = copies.get(x['sequence_key'], 0) + 1
    for x in rows:
        x['key_copies'] = copies[x['sequence_key']]
        x['overlapping_candidates'] = []
    for i, x in enumerate(rows):
        for y in rows[i + 1:]:
            if y['contig'] != x['contig'] or y['start'] > x['end']:
                break
            x['overlapping_candidates'].append(y['inverton'])
            y['overlapping_candidates'].append(x['inverton'])
    for x in rows:
        x['overlapping_candidates'] = '; '.join(x['overlapping_candidates'])
        del x['_order']
    return rows


def write_tsv(path, rows, cols):
    tmp = path + '.tmp'
    with open(tmp, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    os.replace(tmp, path)


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        return [dict(zip(hdr, line.rstrip('\n').split('\t'))) for line in fh if line.strip()]


def median(xs):
    xs = sorted(xs)
    if not xs:
        return None
    k = len(xs) // 2
    return xs[k] if len(xs) % 2 else (xs[k - 1] + xs[k]) / 2


def cmd_annotate(a):
    cons = read_fasta(a.consensus)
    feats = load_features(a.gff)
    ann = annotation_seqs(a.gff, a.fna) if feats else []
    c2a = contig_map(cons, ann) if feats else {}
    ratio = read_ratio(a.ratio)
    geno = None
    if a.genotype:
        geno = {r['inverton']: r for r in read_tsv(a.genotype)}
        need = {r['Inverton'] for r in ratio if int(r['ReverseReadCount']) >= 1}
        if set(geno) != need:
            raise SystemExit(f'{a.genotype} does not cover exactly the candidates with PhaVa inverted reads')
    rows = annotate_rows(a.sample, ratio, cons, feats, c2a, a.min_reverse_reads,
                         a.min_reverse_fraction, a.promoter_window, geno)
    write_tsv(a.out, rows, COLUMNS)
    sub = json.load(open(a.subsample_json)) if a.subsample_json and os.path.exists(a.subsample_json) else {}
    called = [r for r in rows if r['status'] == 'called']
    med = median([r['counted_reads'] for r in rows])
    phava_called = sum(1 for r in rows if r['phava_reverse_reads'] >= a.min_reverse_reads and
                       r['phava_reverse_reads'] / (r['phava_forward_reads'] + r['phava_reverse_reads']) >= a.min_reverse_fraction)
    summ = {'sample': a.sample, 'candidates': len(rows), 'spanned': sum(r['counted_reads'] > 0 for r in rows),
            'with_reverse_reads': sum(r['reverse_reads'] > 0 for r in rows),
            'regenotyped': sum(r['counts_from'] == 'regenotyped' for r in rows), 'phava_called': phava_called,
            'called': len(called),
            'called_intergenic': sum(r['class'] == 'intergenic' for r in called),
            'called_gene_boundary': sum(r['class'] == 'gene_boundary' for r in called),
            'called_intragenic': sum(r['class'] == 'intragenic' for r in called),
            'called_with_promoter_target': sum(bool(r['promoter_targets']) and r['promoter_targets'] != 'NA'
                                               for r in called),
            'median_counted_reads': 'NA' if med is None else (int(med) if med == int(med) else med),
            'reads_used': sub.get('used_reads', 'NA'),
            'depth_used': f"{sub['used_depth']:.1f}" if 'used_depth' in sub else 'NA',
            'subsampled': ('yes' if sub.get('subsampled') else 'no') if sub else 'NA',
            'min_reverse_reads': a.min_reverse_reads, 'min_reverse_fraction': a.min_reverse_fraction,
            'annotation': os.path.basename(a.gff) if feats and c2a else 'none',
            'phava_version': a.phava_version, 'tools_version': VERSION}
    write_tsv(a.summary, [summ], SUMMARY_COLUMNS)
    print(json.dumps({k: summ[k] for k in ('candidates', 'spanned', 'with_reverse_reads', 'called')}))


# ------------------------------------------------------------------------ collect
def cmd_collect(a):
    samples = [s for s in a.samples.split(',') if s]
    per, summ = {}, []
    for s in samples:
        inv = os.path.join(a.dir, s, 'invertons.tsv')
        sm = os.path.join(a.dir, s, 'summary.tsv')
        if not (os.path.exists(inv) and os.path.exists(sm)):
            print(f'warning: no PhaVa tables for {s}', file=sys.stderr)
            continue
        per[s] = read_tsv(inv)
        summ += read_tsv(sm)
    by_key = {}
    for s, rows in per.items():
        for r in rows:
            by_key.setdefault(r['sequence_key'], []).append(r)
    called = []
    for s, rows in per.items():
        for r in rows:
            if r['status'] != 'called':
                continue
            others = [o for o in by_key[r['sequence_key']] if o['sample'] != s]
            r = dict(r)
            r['in_other_samples'] = '; '.join(
                f"{o['sample']} {o['status']} {o['reverse_fraction']} ({o['reverse_reads']}/{o['counted_reads']})"
                for o in sorted(others, key=lambda o: (samples.index(o['sample']), int(o['start']))))
            called.append(r)
    write_tsv(os.path.join(a.dir, 'summary.tsv'), summ, SUMMARY_COLUMNS)
    write_tsv(os.path.join(a.dir, 'invertons_called.tsv'), called, COLUMNS + COMBINED_EXTRA)
    write_tsv(os.path.join(a.dir, 'invertons_all.tsv'), [r for s in per for r in per[s]], COLUMNS)
    print(json.dumps({'samples': len(per), 'called': len(called)}))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', action='version', version=VERSION)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('gff')
    p.add_argument('--gff', required=True)
    p.add_argument('--fna', default='')
    p.add_argument('--consensus', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--json', default='')
    p.add_argument('--transl-table', type=int, default=11)
    p.set_defaults(func=cmd_gff)
    p = sp.add_parser('subsample')
    p.add_argument('--fastq', required=True)
    p.add_argument('--genome', required=True)
    p.add_argument('--max-depth', type=float, required=True)
    p.add_argument('--seed', default='phava1')
    p.add_argument('--out', required=True)
    p.add_argument('--json', required=True)
    p.set_defaults(func=cmd_subsample)
    p = sp.add_parser('genotype')
    p.add_argument('--ratio', required=True)
    p.add_argument('--consensus', required=True)
    p.add_argument('--bam', required=True, help='indexed mapping of the reads to the consensus (s13)')
    p.add_argument('--fastq', required=True, help='the reads PhaVa used')
    p.add_argument('--subsample-json', default='')
    p.add_argument('--flank', type=int, default=1000)
    p.add_argument('--max-mismatch', type=float, default=0.15)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--tmp', default='')
    p.add_argument('--out', required=True)
    p.set_defaults(func=cmd_genotype)
    p = sp.add_parser('annotate')
    p.add_argument('--ratio', required=True)
    p.add_argument('--consensus', required=True)
    p.add_argument('--gff', default='')
    p.add_argument('--fna', default='')
    p.add_argument('--sample', required=True)
    p.add_argument('--min-reverse-reads', type=int, default=3)
    p.add_argument('--min-reverse-fraction', type=float, default=0.03)
    p.add_argument('--promoter-window', type=int, default=300)
    p.add_argument('--subsample-json', default='')
    p.add_argument('--genotype', default='', help='genotype.tsv: re-counted candidates (default: PhaVa counts)')
    p.add_argument('--phava-version', default='')
    p.add_argument('--out', required=True)
    p.add_argument('--summary', required=True)
    p.set_defaults(func=cmd_annotate)
    p = sp.add_parser('collect')
    p.add_argument('--dir', required=True)
    p.add_argument('--samples', required=True)
    p.set_defaults(func=cmd_collect)
    a = ap.parse_args()
    a.func(a)


if __name__ == '__main__':
    main()
