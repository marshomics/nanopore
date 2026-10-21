#!/usr/bin/env python3
"""Summarise the s18b_sv validation: one JSON per simulated condition -> tables + figure."""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RES, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(RES, "*.json")))]
df = pd.DataFrame(rows)
df["depth"] = df["depth"].astype(int)
df["cons"] = np.where(df["consensus"].str.contains("deletion"), "B", "A")
df["panel"] = df["pair"] + np.where(df["pair"] == "DEL", " (consensus " + df["cons"] + ")", "")
f = df["true_alt_fraction"]
df["mix"] = (df["nominal_alt_fraction"] > 0) & (df["nominal_alt_fraction"] < 1)


def called_mixture(st):
    return st == "mixture"


df["known_called_mixture"] = df["known_status"].map(called_mixture)
df["known_j_called_mixture"] = df["known_j_status"].map(called_mixture)
df["disc_reported"] = df["disc_found"].fillna(False) & df["disc_status"].isin(["mixture", "consensus_discrepancy"])
keep = ["cond", "panel", "depth", "nominal_alt_fraction", "true_alt_fraction", "known_status", "known_spanning",
        "known_alt_reads", "known_alt_fraction", "known_ci_lo", "known_ci_hi", "known_abs_error", "known_ci_covers_truth",
        "known_reads_checked", "known_reads_misassigned", "known_j_status", "known_j_fraction", "known_j_ci_lo",
        "known_j_ci_hi", "known_j_abs_error", "known_j_ci_covers_truth", "disc_found", "disc_kind", "disc_size",
        "disc_status", "disc_method", "disc_alt_fraction", "disc_truth", "disc_abs_error", "disc_j_fraction",
        "disc_j_abs_error", "disc_candidates", "disc_false_reported", "disc_false_list", "elapsed_s", "disc_elapsed_s"]
df.sort_values(["panel", "depth", "nominal_alt_fraction"])[[c for c in keep if c in df.columns]].to_csv(
    os.path.join(OUT, "validation_results.tsv"), sep="\t", index=False)

# ---- summary per panel x depth
summ = []
for (pn, d), g in df.groupby(["panel", "depth"]):
    m = g[g["mix"]]
    mm = m[(m["nominal_alt_fraction"] >= 0.05) & (m["nominal_alt_fraction"] <= 0.95)]
    pure = g[~g["mix"]]
    det = m[m["known_called_mixture"]]
    summ.append({
        "event": pn, "depth": d, "conditions": len(g),
        "median_spanning_reads": int(g["known_spanning"].median()),
        "spanning_MAE_pp_5to95": round(100 * mm["known_abs_error"].mean(), 2),
        "spanning_bias_pp_5to95": round(100 * (mm["known_alt_fraction"] - mm["true_alt_fraction"]).mean(), 2),
        "spanning_max_err_pp_5to95": round(100 * mm["known_abs_error"].max(), 2),
        "spanning_CI_covers_truth": f"{int(m['known_ci_covers_truth'].sum())}/{len(m)}",
        "junction_MAE_pp_5to95": round(100 * mm["known_j_abs_error"].mean(), 2),
        "junction_bias_pp_5to95": round(100 * (mm["known_j_fraction"] - mm["true_alt_fraction"]).mean(), 2),
        "junction_CI_covers_truth": f"{int(m['known_j_ci_covers_truth'].sum())}/{len(m)}",
        "lowest_fraction_called_mixture": (f"{100 * det['nominal_alt_fraction'].min():g}%" if len(det) else "none"),
        "mixtures_called": f"{len(det)}/{len(m)}",
        "false_mixture_calls_in_pure": int(pure["known_called_mixture"].sum()),
        "reads_misassigned": f"{int(g['known_reads_misassigned'].sum())}/{int(g['known_reads_checked'].sum())}",
        "discovery_reported_5to95": f"{int(mm['disc_reported'].sum())}/{len(mm)}",
        "discovery_MAE_pp_5to95": round(100 * mm.loc[mm['disc_reported'], 'disc_abs_error'].mean(), 2)
        if mm["disc_reported"].any() else None,
        "false_discovered_calls": int(g["disc_false_reported"].sum()),
        "median_runtime_s": int((g["elapsed_s"] + g["disc_elapsed_s"].fillna(0)).median()),
    })
S = pd.DataFrame(summ)
mixall = df[df["mix"]]
pooled = {"conditions": len(df), "mixtures": len(mixall),
          "spanning_CI_covers_truth": f"{int(mixall['known_ci_covers_truth'].sum())}/{len(mixall)}",
          "junction_CI_covers_truth": f"{int(mixall['known_j_ci_covers_truth'].sum())}/{len(mixall)}",
          "false_mixture_calls_in_pure": int(df.loc[~df['mix'], 'known_called_mixture'].sum()),
          "pure_conditions": int((~df['mix']).sum()),
          "reads_misassigned": f"{int(df['known_reads_misassigned'].sum())}/{int(df['known_reads_checked'].sum())}",
          "false_discovered_calls": int(df['disc_false_reported'].sum())}
json.dump(pooled, open(os.path.join(OUT, "validation_pooled.json"), "w"), indent=1)
print(json.dumps(pooled))
S.to_csv(os.path.join(OUT, "validation_summary.tsv"), sep="\t", index=False)
print(S.to_string(index=False))

# ---- figure: estimated vs true fraction, spanning reads; zoom on minor fractions; junction method
plt.rcParams.update({"font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "legend.fontsize": 6, "axes.spines.top": False,
                     "axes.spines.right": False, "font.family": "DejaVu Sans"})
panels = ["DEL (consensus A)", "DEL (consensus B)", "INV", "DUP"]
titles = {"DEL (consensus A)": "6.5-kb deletion;\nfull version assembled",
          "DEL (consensus B)": "6.5-kb deletion;\ndeletion version assembled",
          "INV": "3-kb inversion", "DUP": "5-kb tandem\nduplication"}
depths = sorted(df["depth"].unique())
cmap = plt.get_cmap("Blues")
col = {d: cmap(0.45 + 0.5 * i / max(1, len(depths) - 1)) for i, d in enumerate(depths)}
fig, axes = plt.subplots(3, 4, figsize=(7.2, 5.6), sharex="row", sharey="row")
for j, pn in enumerate(panels):
    g = df[df["panel"] == pn]
    for row, (est, lo, hi, lim) in enumerate([("known_alt_fraction", "known_ci_lo", "known_ci_hi", 1.0),
                                               ("known_alt_fraction", "known_ci_lo", "known_ci_hi", 0.12),
                                               ("known_j_fraction", "known_j_ci_lo", "known_j_ci_hi", 1.0)]):
        ax = axes[row, j]
        ax.plot([0, lim], [0, lim], color="0.6", lw=0.6, zorder=1)
        if row == 1:
            ax.axhline(0.01, color="0.4", lw=0.6, ls=":", zorder=1)
        for i, d in enumerate(depths):
            h = g[(g["depth"] == d) & g[est].notna()]
            if row == 1:
                h = h[h["true_alt_fraction"] <= lim]
            if h.empty:
                continue
            x = h["true_alt_fraction"] + (i - 1) * lim * 0.012
            y = h[est].clip(upper=lim * 1.15)
            yerr = np.vstack([(h[est] - h[lo]).clip(lower=0), (h[hi].clip(upper=lim * 1.15) - h[est]).clip(lower=0)])
            ax.errorbar(x, y, yerr=yerr, fmt="o", ms=2.6, lw=0.7, color=col[d], mec=col[d], capsize=0,
                        label=f"{d}×", zorder=3)
        ax.set_xlim(-0.03 * lim, lim * 1.04)
        ax.set_ylim(-0.03 * lim, lim * 1.08)
        if row == 0:
            ax.set_title(titles[pn], loc="left")
        if j == 0:
            ax.set_ylabel(["Estimated fraction\n(spanning reads)", "Estimated fraction,\nminor version (spanning)",
                           "Estimated fraction\n(breakpoint junctions)"][row])
        if row == 2 or row == 1:
            ax.set_xlabel("True fraction" if row == 2 else "True fraction (0–12%)")
        if row == 1:
            tk = [0, 0.04, 0.08, 0.12]
            ax.set_xticks(tk)
            ax.set_yticks(tk)
            ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
            ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
        else:
            ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
axes[0, 0].legend(title="Read depth", frameon=False, loc="upper left", handletextpad=0.2, borderaxespad=0.1)
axes[1, 0].text(0.118, 0.0112, "1% floor", ha="right", va="bottom", fontsize=6, color="0.3")
for ax, letter in zip(axes[:, 0], "abc"):
    ax.text(-0.42, 1.08, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
fig.tight_layout(h_pad=1.0, w_pad=0.5)
fig.savefig(os.path.join(OUT, "validation_accuracy.png"), dpi=300)
fig.savefig(os.path.join(OUT, "validation_accuracy.pdf"))
print("figure written")
