#!/usr/bin/env python3
"""Describe (1) the parent reads found in only one of S3 / S1 / P and (2) the size of the ML/MM
differences between S3 and S1.  Usage: inspect_diffs.py S3.bam S1.bam P.bam F1,F2,...  -> JSON"""
import collections, json, sys
import pysam
s3p, s1p, pp, fns = sys.argv[1], sys.argv[2], sys.argv[3], set(sys.argv[4].split(","))

def load(path):
    recs = {}
    with pysam.AlignmentFile(path, check_sq=False) as f:
        for r in f.fetch(until_eof=True):
            if r.has_tag("fn") and r.get_tag("fn") not in fns:
                continue
            recs[r.query_name] = r
    return recs
S3, S1, P = load(s3p), load(s1p), load(pp)
def parent(r):
    return r.get_tag("pi") if r.has_tag("pi") else r.query_name
def by_parent(recs):
    d = collections.defaultdict(list)
    for r in recs.values():
        d[parent(r)].append(r)
    return d
b3, b1, bp = by_parent(S3), by_parent(S1), by_parent(P)
def describe(rs):
    return [{"read": r.query_name, "len": r.query_length, "qs": round(float(r.get_tag("qs")), 2) if r.has_tag("qs") else None,
             "ns": r.get_tag("ns") if r.has_tag("ns") else None, "ts": r.get_tag("ts") if r.has_tag("ts") else None,
             "has_pi": r.has_tag("pi"), "fn": r.get_tag("fn") if r.has_tag("fn") else None} for r in rs]
out = {"parents_only_S3_vs_P": {}, "parents_only_P_vs_S3": {}}
for p in sorted(set(b3) - set(bp)):
    out["parents_only_S3_vs_P"][p] = {"S3": describe(b3[p]), "in_S1": p in b1, "S1": describe(b1.get(p, []))}
for p in sorted(set(bp) - set(b3)):
    out["parents_only_P_vs_S3"][p] = {"P": describe(bp[p]), "in_S1": p in b1}
# size of the ML differences between S3 and S1 (same reads)
mags, npos, mm_diff = [], [], 0
for k in S3:
    a, b = S3[k], S1.get(k)
    if b is None:
        continue
    if a.get_tag("MM") != b.get_tag("MM"):
        mm_diff += 1
    la, lb = list(a.get_tag("ML")), list(b.get_tag("ML"))
    if la != lb:
        if len(la) == len(lb):
            d = [abs(x - y) for x, y in zip(la, lb) if x != y]
            mags.append(max(d)); npos.append(len(d))
        else:
            npos.append(-1)
out["ml_differing_reads"] = len(npos)
out["ml_same_length_reads"] = sum(1 for x in npos if x >= 0)
out["ml_positions_differing_per_read"] = sorted(x for x in npos if x >= 0)
out["ml_max_abs_difference_per_read_of_255"] = sorted(mags)
out["mm_differing_reads"] = mm_diff
json.dump(out, sys.stdout, indent=1, sort_keys=True)
print()
