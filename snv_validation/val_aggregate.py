#!/usr/bin/env python3
"""val_aggregate.py - collect the three small-variant validations into evidence folders.

  polval  rule R1 tables: summary.tsv (Dref, H), pairs.tsv (Drep), injections.tsv, raw/
  snvval  rule R2 tables: summary.tsv (vcfdist counts, pooled P/R/F1), strata.tsv, raw/
  minval  rule R3 tables: conditions.tsv, summary.tsv, sensitivity.tsv, raw/

Metrics here are computed independently of the gate-check scripts, which recompute them
from the raw files and must agree.
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys

S = '/ebio/abt3_scratch/jmarsh/software/science/snv_build'
SAMPLES = ['Cminuta01', 'Cminuta02', 'Tube14_T002_conj', 'Tube15_T002_conj']
PAIRS = {'Cminuta': ('Cminuta01', 'Cminuta02'), 'T002': ('Tube14_T002_conj', 'Tube15_T002_conj')}
STATES = ['draft', 'd50', 'd100', 'd200', 'd400']


def tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith('##')), delimiter='\t'))


def write(path, cols, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')


def state_fasta(sample, state):
    d = f'{S}/polval/{sample}'
    return f'{d}/draft.fasta' if state == 'draft' else f'{d}/polished_{state}.fasta'


# ---------------------------------------------------------------------- polval
def polval(out):
    raw = f'{out}/raw'
    for sub in ('dref', 'drep', 'holdout', 'inject', 'changes'):
        os.makedirs(f'{raw}/{sub}', exist_ok=True)
    rows = []
    for s in SAMPLES:
        for x in STATES:
            src = f'{S}/polval/{s}/dref_{x}.tsv'
            shutil.copyfile(src, f'{raw}/dref/{s}.{x}.tsv')
            d = tsv(src)
            small = [r for r in d if int(r['length']) <= 50]
            vcf = f'{S}/polval/{s}/c3_ho_{x}.vcf'
            shutil.copyfile(vcf, f'{raw}/holdout/{s}.{x}.vcf')
            npass = sum(1 for l in open(vcf) if not l.startswith('#') and l.split('\t')[6] == 'PASS')
            rows.append(dict(sample=s, state=x, Dref=len(small), Dref_snp=sum(r['kind'] == 'SNP' for r in small),
                             Dref_indel=sum(r['kind'] in ('INS', 'DEL') for r in small),
                             Dref_large=len(d) - len(small), H=npass))
    write(f'{out}/summary.tsv', ['sample', 'state', 'Dref', 'Dref_snp', 'Dref_indel', 'Dref_large', 'H'], rows)
    prow = []
    for p, (a, b) in PAIRS.items():
        for x in STATES:
            dst = f'{raw}/drep/{p}.{x}.tsv'
            subprocess.run([sys.executable, f'{S}/bin/polish_tools.py', 'diff', '--query', state_fasta(b, x),
                            '--ref', state_fasta(a, x), '--out', dst, '--threads', '8'], check=True,
                           stdout=subprocess.DEVNULL)
            os.remove(dst + '.json')
            d = tsv(dst)
            small = [r for r in d if int(r['length']) <= 50]
            prow.append(dict(pair=p, state=x, Drep=len(small), Drep_snp=sum(r['kind'] == 'SNP' for r in small),
                             Drep_indel=sum(r['kind'] in ('INS', 'DEL') for r in small), Drep_large=len(d) - len(small)))
    write(f'{out}/pairs.tsv', ['pair', 'state', 'Drep', 'Drep_snp', 'Drep_indel', 'Drep_large'], prow)
    irow = []
    for s in SAMPLES:
        for dp in (50, 100, 200, 400):
            a_ = f'{S}/polval/{s}/inj_d{dp}.assess.tsv'
            n_ = f'{S}/polval/{s}/inj_d{dp}.new.tsv'
            shutil.copyfile(a_, f'{raw}/inject/{s}.d{dp}.tsv')
            shutil.copyfile(n_, f'{raw}/inject/{s}.d{dp}.new.tsv')
            shutil.copyfile(f'{S}/polval/{s}/changes_d{dp}.tsv', f'{raw}/changes/{s}.d{dp}.tsv')
            t = tsv(a_)
            irow.append(dict(sample=s, depth=dp, injected_snp=sum(r['kind'] == 'SNP' for r in t),
                             injected_indel=sum(r['kind'] == 'HP' for r in t),
                             restored_snp=sum(int(r['restored']) for r in t if r['kind'] == 'SNP'),
                             restored_indel=sum(int(r['restored']) for r in t if r['kind'] == 'HP'),
                             new_diffs=len(tsv(n_))))
    write(f'{out}/injections.tsv', ['sample', 'depth', 'injected_snp', 'injected_indel', 'restored_snp',
                                    'restored_indel', 'new_diffs'], irow)
    for s in SAMPLES:
        shutil.copyfile(f'{S}/polval/{s}/sub.subsample.json', f'{raw}/{s}.subsample.json')
    print(json.dumps({'summary': len(rows), 'pairs': len(prow), 'injections': len(irow)}))


# ---------------------------------------------------------------------- snvval
def load_bed(path):
    out = {}
    for l in open(path):
        f = l.split()
        if len(f) >= 3:
            out.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return out


def in_bed(bed, c, p0):
    return any(s <= p0 < e for s, e in bed.get(c, []))


def snvval(out):
    raw = f'{out}/raw'
    os.makedirs(raw, exist_ok=True)
    rows, strata = [], []
    for s in SAMPLES:
        W = f'{S}/snvval/{s}'
        kind = {}
        for l in open(f'{W}/truth.vcf'):
            if l.startswith('#'):
                continue
            f = l.rstrip('\n').split('\t')
            info = dict(kv.split('=', 1) for kv in f[7].split(';') if '=' in kv)
            kind[(f[0], int(f[1]))] = info.get('KIND', '?')
        rep = load_bed(f'{W}/R.repeats.bed')
        for dp in (10, 25, 50, 100):
            for cf in ('C3P', 'C3AF', 'DP', 'INT', 'UNI'):
                src = f'{W}/vd_{cf}_d{dp}/precision-recall-summary.tsv'
                shutil.copyfile(src, f'{raw}/{s}.{dp}.{cf}.precision-recall-summary.tsv')
                v = {r['VAR_TYPE']: r for r in tsv(src) if r['THRESHOLD'] == 'NONE'}
                c = {}
                for vt, pre in (('SNP', 'snp'), ('INDEL', 'indel')):
                    c[f'{pre}_ttp'] = int(v[vt]['TRUTH_TP'])
                    c[f'{pre}_qtp'] = int(v[vt]['QUERY_TP'])
                    c[f'{pre}_fp'] = int(v[vt]['QUERY_FP'])
                    c[f'{pre}_fn'] = int(v[vt]['TRUTH_FN'])
                ttp, qtp = c['snp_ttp'] + c['indel_ttp'], c['snp_qtp'] + c['indel_qtp']
                fp, fn = c['snp_fp'] + c['indel_fp'], c['snp_fn'] + c['indel_fn']
                prec = qtp / (qtp + fp) if qtp + fp else 1.0
                rec = ttp / (ttp + fn) if ttp + fn else 1.0
                f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
                rows.append(dict(sample=s, depth=dp, config=cf, **c, precision=f'{prec:.6f}', recall=f'{rec:.6f}',
                                 f1=f'{f1:.6f}'))
                # strata from per-variant truth decisions (POS in truth.tsv is 0-based)
                tt = tsv(f'{W}/vd_{cf}_d{dp}/truth.tsv')
                agg = {}
                for r in tt:
                    k = kind.get((r['CONTIG'], int(r['POS']) + 1), 'other')
                    rp = 'repeat' if in_bed(rep, r['CONTIG'], int(r['POS'])) else 'unique'
                    a = agg.setdefault((k, rp), [0, 0])
                    a[1] += 1
                    a[0] += r['ERRTYPE'].strip() == 'TP'
                qq = tsv(f'{W}/vd_{cf}_d{dp}/query.tsv')
                fps = {}
                for r in qq:
                    if r['ERRTYPE'].strip() == 'FP':
                        rp = 'repeat' if in_bed(rep, r['CONTIG'], int(r['POS'])) else 'unique'
                        fps[rp] = fps.get(rp, 0) + 1
                for (k, rp), (tp, n) in sorted(agg.items()):
                    strata.append(dict(sample=s, depth=dp, config=cf, kind=k, region=rp, truth=n, tp=tp,
                                       recall=f'{tp / n:.5f}' if n else '', fp_in_region=fps.get(rp, 0)))
    write(f'{out}/summary.tsv', ['sample', 'depth', 'config', 'snp_ttp', 'snp_qtp', 'snp_fp', 'snp_fn', 'indel_ttp',
                                 'indel_qtp', 'indel_fp', 'indel_fn', 'precision', 'recall', 'f1'], rows)
    write(f'{out}/strata.tsv', ['sample', 'depth', 'config', 'kind', 'region', 'truth', 'tp', 'recall', 'fp_in_region'],
          strata)
    print(json.dumps({'summary': len(rows), 'strata': len(strata)}))


# ---------------------------------------------------------------------- minval
MSAMPLES = ['Cminuta01', 'Tube14_T002_conj']
DEPTHS = [100, 500, 2000]
FRACS = ['0.0025', '0.005', '0.01', '0.02', '0.05', '0.10', '0.20', '0.40']


def minval(out):
    raw = f'{out}/raw'
    os.makedirs(raw, exist_ok=True)
    conds, summ, sens = [], [], []
    for s in MSAMPLES:
        W = f'{S}/minval/{s}'
        for D in DEPTHS:
            ctrl = f'{s}.d{D}.ctrl'
            shutil.copyfile(f'{W}/out/ctrl_d{D}/sites.tsv', f'{raw}/{ctrl}.calls.tsv')
            shutil.copyfile(f'{W}/out/ctrl_d{D}/summary.tsv', f'{raw}/{ctrl}.summary.tsv')
            conds.append(dict(cond=ctrl, sample=s, depth=D, fraction=0.0))
            cpos = {(r['contig'], int(r['pos'])) for r in tsv(f'{W}/out/ctrl_d{D}/sites.tsv')}
            for f in FRACS:
                cond = f'{s}.d{D}.f{f}'
                o = f'{W}/out/f{f}_d{D}'
                shutil.copyfile(f'{o}/sites.tsv', f'{raw}/{cond}.calls.tsv')
                shutil.copyfile(f'{o}.truth.tsv', f'{raw}/{cond}.truth.tsv')
                shutil.copyfile(f'{o}/summary.tsv', f'{raw}/{cond}.summary.tsv')
                if os.path.exists(f'{o}/groups.tsv'):
                    shutil.copyfile(f'{o}/groups.tsv', f'{raw}/{cond}.groups.tsv')
                conds.append(dict(cond=cond, sample=s, depth=D, fraction=float(f)))
                truth = tsv(f'{o}.truth.tsv')
                calls = {(r['contig'], int(r['pos']), r['alt']): r for r in tsv(f'{o}/sites.tsv')}
                tpos = {(t['contig'], int(t['pos'])) for t in truth}
                n_sp = n_hit = n_cov = n_ok = masked_called = 0
                errs, by_set = [], {}
                for t in truth:
                    k = (t['contig'], int(t['pos']), t['alt'])
                    hit = calls.get(k)
                    st = t['set']
                    b = by_set.setdefault(st, {'n': 0, 'called': 0, 'linked': 0})
                    b['n'] += 1
                    if hit is not None:
                        b['called' if hit['status'] in ('called', 'majority_alt') else 'linked'] += 1
                    if st == 'masked':
                        masked_called += hit is not None
                        continue
                    n_sp += 1
                    if hit is None:
                        continue
                    n_hit += 1
                    tf, ef = float(t['truth_fraction']), float(hit['fraction'])
                    errs.append(ef - tf)
                    if tf >= 0.01:
                        n_cov += 1
                        n_ok += float(hit['ci_low']) <= tf <= float(hit['ci_high'])
                false = background = 0
                for (c_, p_, a_), r in calls.items():
                    if (c_, p_) in tpos:
                        continue
                    if (c_, p_) in cpos:
                        background += 1
                    else:
                        false += 1
                summ.append(dict(cond=cond, n_spiked=n_sp, n_called_spiked=n_hit, n_calls=len(calls),
                                 n_background=background, n_false=false, n_cov_eligible=n_cov, n_cov_ok=n_ok,
                                 masked_called=masked_called,
                                 mae=f'{(sum(abs(e) for e in errs) / len(errs)) if errs else 0.0:.8f}',
                                 bias=f'{(sum(errs) / len(errs)) if errs else 0.0:.8f}'))
                med_truth = sorted(float(t['truth_fraction']) for t in truth if t['set'] != 'masked')
                for st, b in by_set.items():
                    sens.append(dict(cond=cond, sample=s, depth=D, fraction=f, set=st, n=b['n'],
                                     called=b['called'], linked=b['linked'],
                                     median_truth_fraction=f'{med_truth[len(med_truth) // 2]:.5f}'))
    write(f'{out}/conditions.tsv', ['cond', 'sample', 'depth', 'fraction'], conds)
    write(f'{out}/summary.tsv', ['cond', 'n_spiked', 'n_called_spiked', 'n_calls', 'n_background', 'n_false',
                                 'n_cov_eligible', 'n_cov_ok', 'masked_called', 'mae', 'bias'], summ)
    write(f'{out}/sensitivity.tsv', ['cond', 'sample', 'depth', 'fraction', 'set', 'n', 'called', 'linked',
                                     'median_truth_fraction'], sens)
    print(json.dumps({'conditions': len(conds), 'summary': len(summ)}))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('what', choices=['polval', 'snvval', 'minval'])
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    {'polval': polval, 'snvval': snvval, 'minval': minval}[a.what](a.out)
