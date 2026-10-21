#!/usr/bin/env python3
"""Which reads differ between S3 and S1, between S1' and S1 and between S3 and S1' (records reduced as
in compare_bams.py), and how far these sets overlap; for each pair, which fields differ.
Descriptive only (not an R5 criterion).  Usage: overlap_diffs.py S3.bam S1.bam S1b.bam -> JSON"""
import collections, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_bams import load

def diff(x, y):
    out, fields = {}, collections.Counter()
    for k in set(x) & set(y):
        a, b = x[k], y[k]
        d = [f for f in ("flag", "seq", "qual") if a[f] != b[f]]
        d += ["tag:" + t for t in sorted(set(a["tags"]) | set(b["tags"])) if a["tags"].get(t) != b["tags"].get(t)]
        if d:
            out[k] = d; fields.update(d)
    return out, dict(fields)

S3, S1, S1b = (load(p)[1] for p in sys.argv[1:4])
d31, f31 = diff(S3, S1)
db1, fb1 = diff(S1b, S1)
d3b, f3b = diff(S3, S1b)
a, b, c = set(d31), set(db1), set(d3b)
out = {
    "reads": {"S3": len(S3), "S1": len(S1), "S1b": len(S1b)},
    "same_read_ids": set(S3) == set(S1) == set(S1b),
    "differing": {"S3_vs_S1": len(a), "S1b_vs_S1": len(b), "S3_vs_S1b": len(c)},
    "fields": {"S3_vs_S1": f31, "S1b_vs_S1": fb1, "S3_vs_S1b": f3b},
    "overlap_S3vsS1_and_S1bvsS1": len(a & b),
    "S3vsS1_only": len(a - b),
    "S1bvsS1_only": len(b - a),
    "differing_in_any_pair": len(a | b | c),
    "S3_vs_S1_reads": sorted(a), "S1b_vs_S1_reads": sorted(b), "S3_vs_S1b_reads": sorted(c),
}
json.dump(out, sys.stdout, indent=1, sort_keys=True)
print()
