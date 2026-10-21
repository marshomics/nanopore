#!/usr/bin/env python3
"""Rule R4 metrics from the pulled validation tables: metrics.tsv (one row per sample and
fraction) and val_invertons.tsv (one row per spiked inverton and fraction); the figure is
plot_accuracy.py. argv: val-dir"""
import os
import statistics
import sys

import pandas as pd

V = sys.argv[1]
SAMPLES = ['Cminuta01', 'Tube14_T002_conj']
FR = ['0.05', '0.2', '0.5']


def rd(p):
    return pd.read_csv(p, sep='\t', dtype=str, keep_default_na=False)


met, per = [], []
for s in SAMPLES:
    d = os.path.join(V, s)
    spiked = rd(os.path.join(d, 'spiked.tsv'))
    sp = set(spiked['inverton'])
    ctrl = rd(os.path.join(d, 'control', 'invertons.tsv')).set_index('inverton')
    ctrl_called = set(ctrl.index[ctrl['status'] == 'called']) - sp
    truth = rd(os.path.join(d, 'truth.tsv'))
    truth['fraction'] = truth['fraction'].astype(float).map(lambda x: f'{x:g}')
    for f in FR:
        cond = rd(os.path.join(d, f'f{f}', 'invertons.tsv')).set_index('inverton')
        ratio = rd(os.path.join(d, f'f{f}', 'phava_ratio.tsv')).set_index('Inverton')
        t = truth[truth['fraction'] == f].set_index('inverton')
        diffs, ncalled, pdiffs, pcalled = [], 0, [], 0
        for inv in spiked['inverton']:
            c = cond.loc[inv]
            F, R = int(c['forward_reads']), int(c['reverse_reads'])
            P = R / (F + R) if F + R else float('nan')
            fr = float(t.loc[inv, 'f_real'])
            called = c['status'] == 'called'
            if called:
                ncalled += 1
                diffs.append(abs(P - fr))
            pF, pR = int(ratio.loc[inv, 'ForwardReadCount']), int(ratio.loc[inv, 'ReverseReadCount'])
            pP = float(ratio.loc[inv, 'ProportionReverse'])
            if pR >= 3 and pR / (pF + pR) >= 0.03:
                pcalled += 1
                pdiffs.append(abs(pP - fr))
            per.append({'sample': s, 'fraction': f, 'inverton': inv, 'segment_bp': c['segment_bp'], 'class': c['class'],
                        'spanning_reads': t.loc[inv, 'spanning_reads'], 'edited_reads': t.loc[inv, 'edited_reads'],
                        'f_real': f'{fr:.4f}', 'forward_reads': F, 'reverse_reads': R, 'reverse_fraction': f'{P:.4f}',
                        'ci_low': c['ci_low'], 'ci_high': c['ci_high'], 'status': c['status'],
                        'abs_diff': f'{abs(P - fr):.4f}', 'phava_forward_reads': pF, 'phava_reverse_reads': pR,
                        'phava_proportion_reverse': f'{pP:.4f}', 'control_forward_reads': ctrl.loc[inv, 'forward_reads']})
        unsp = set(cond.index[cond['status'] == 'called']) - sp
        met.append({'sample': s, 'fraction': f, 'spiked': '20', 'called_spiked': str(ncalled),
                    'within_0.05': str(sum(x <= 0.05 for x in diffs)),
                    'median_abs_diff': f'{statistics.median(diffs):.4f}' if diffs else 'NA',
                    'max_abs_diff': f'{max(diffs):.4f}' if diffs else 'NA',
                    'phava_called_spiked': str(pcalled), 'phava_within_0.05': str(sum(x <= 0.05 for x in pdiffs)),
                    'phava_median_abs_diff': f'{statistics.median(pdiffs):.4f}' if pdiffs else 'NA',
                    'called_unspiked': str(len(unsp)), 'control_called': str(len(ctrl_called)),
                    'spill_over_equal': 'yes' if unsp == ctrl_called else 'no',
                    'control_clean': 'yes' if all(int(ctrl.loc[i, 'reverse_reads']) == 0 for i in sp) else 'no'})
pd.DataFrame(met).to_csv(os.path.join(V, 'metrics.tsv'), sep='\t', index=False)
pd.DataFrame(per).to_csv(os.path.join(V, 'val_invertons.tsv'), sep='\t', index=False)
print(pd.DataFrame(met).to_string(index=False))

