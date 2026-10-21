#!/usr/bin/env python3
"""Figure for rule R4: inverted fraction of every spiked inverton against the truth, re-counted
(stage output) and PhaVa's joint count. argv: val_invertons.tsv out.png"""
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False})
x = pd.read_csv(sys.argv[1], sep='\t', dtype=str)
for c in ('f_real', 'reverse_fraction', 'phava_proportion_reverse'):
    x[c] = x[c].astype(float)
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), sharex=True, sharey=True, gridspec_kw={'wspace': 0.08})
names = {'Cminuta01': r'$\mathit{C.\ minuta}$ (Cminuta01)', 'Tube14_T002_conj': r'$\mathit{M.\ smithii}$ (Tube14_T002_conj)'}
for ax, s in zip(axes, ['Cminuta01', 'Tube14_T002_conj']):
    y = x[x['sample'] == s]
    ax.fill_between([0, 0.62], [-0.05, 0.57], [0.05, 0.67], color='0.92', lw=0, zorder=0, label='±0.05 of the truth')
    ax.plot([0, 0.62], [0, 0.62], color='0.35', lw=0.7, zorder=1)
    ax.scatter(y['f_real'], y['phava_proportion_reverse'], s=18, facecolors='none', edgecolors='0.55', lw=0.8, zorder=2,
               label="PhaVa's joint count")
    ax.scatter(y['f_real'], y['reverse_fraction'], s=12, color='#1f5fa8', zorder=3, label='re-counted (stage output)')
    bad = y[(y['phava_proportion_reverse'] - y['f_real']).abs() > 0.05]['inverton'].nunique()
    ax.set_title(names[s], loc='left')
    ax.text(0.98, 0.10, f"PhaVa joint count off by >0.05:\n{bad} of {y['inverton'].nunique()} invertons",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=7, color='0.3')
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 0.92)
    ax.margins(0.04)
    ax.set_xticks([0.05, 0.2, 0.5])
    ax.set_xticklabels(['0.05', '0.2', '0.5'])
    ax.set_xlabel('spiked inverted fraction (edited spanning reads)')
axes[0].set_ylabel('estimated inverted fraction')
axes[0].legend(loc='upper left', frameon=False, handletextpad=0.3)
fig.suptitle("Re-counted fractions match the spiked truth; PhaVa's joint counts overestimate 6 of 40 invertons",
             x=0.01, ha='left', fontsize=8)
fig.subplots_adjust(left=0.08, right=0.99, bottom=0.14, top=0.86)
fig.canvas.draw()
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip() and t.get_visible()]
bad_pairs = [(a.get_text()[:28], b.get_text()[:28]) for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
outside = [t.get_text()[:28] for t, bb in texts if not (fig.bbox.x0 <= bb.x0 and bb.x1 <= fig.bbox.x1 and fig.bbox.y0 <= bb.y0 and bb.y1 <= fig.bbox.y1)]
fig.savefig(sys.argv[2], dpi=300)
print('overlaps', bad_pairs, '| outside', outside)
