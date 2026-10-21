#!/usr/bin/env python3
"""Compare the demultiplexed reads of two pipeline runs, per output file:
s02_demux FASTQs (fastq_raw/) and s02b_demux_bam BAMs (bams/).  Usage: demux_compare.py RUN_A RUN_B"""
import gzip, json, os, sys
import pysam


def fastq_ids(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return [line[1:].split()[0] for i, line in enumerate(fh) if i % 4 == 0]


def bam_ids(path):
    with pysam.AlignmentFile(path, check_sq=False) as f:
        return [r.query_name for r in f.fetch(until_eof=True)]


def collect(d, kind):
    out = {}
    for fn in sorted(os.listdir(d)):
        p = os.path.join(d, fn)
        if not os.path.isfile(p):
            continue
        if kind == "fastq" and fn.endswith((".fastq", ".fastq.gz")):
            out[fn] = fastq_ids(p)
        elif kind == "bam" and fn.endswith(".bam"):
            out[fn] = bam_ids(p)
    return out


def comp(a, b):
    res = {"files_a": len(a), "files_b": len(b),
           "files_only_a": sorted(set(a) - set(b)), "files_only_b": sorted(set(b) - set(a)),
           "reads_a": sum(len(v) for v in a.values()), "reads_b": sum(len(v) for v in b.values()),
           "per_file": {}}
    ok = bool(a) and not res["files_only_a"] and not res["files_only_b"]
    for f in sorted(set(a) & set(b)):
        sa, sb = set(a[f]), set(b[f])
        r = {"reads_a": len(a[f]), "reads_b": len(b[f]), "dup_a": len(a[f]) - len(sa),
             "dup_b": len(b[f]) - len(sb), "only_a": len(sa - sb), "only_b": len(sb - sa)}
        res["per_file"][f] = r
        ok = ok and not (r["only_a"] or r["only_b"] or r["dup_a"] or r["dup_b"])
    res["identical"] = ok
    return res


if __name__ == "__main__":
    ra, rb = sys.argv[1], sys.argv[2]
    out = {"run_a": ra, "run_b": rb,
           "s02_fastq": comp(collect(os.path.join(ra, "fastq_raw"), "fastq"), collect(os.path.join(rb, "fastq_raw"), "fastq")),
           "s02b_bam": comp(collect(os.path.join(ra, "bams"), "bam"), collect(os.path.join(rb, "bams"), "bam"))}
    json.dump(out, sys.stdout, indent=1, sort_keys=True)
    print()
