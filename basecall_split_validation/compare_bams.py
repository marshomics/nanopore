#!/usr/bin/env python3
"""Compare unaligned dorado BAMs read by read (rule R5).

Usage:
  compare_bams.py summary BAM [--fn F1,F2,...]              -> JSON on stdout
  compare_bams.py compare A B [--fn F1,F2,...]              -> JSON on stdout

Each record is reduced to digests of FLAG, SEQ, QUAL and every tag (type and value). With --fn only
records whose fn tag names one of the files are used (for the production BAM, which holds all files).
"""
import argparse, collections, hashlib, json, sys
import pysam


def digest(x):
    return hashlib.sha1(repr(x).encode()).hexdigest()[:16]


def load(path, fns=None):
    recs, dup, seen_total = {}, 0, 0
    with pysam.AlignmentFile(path, check_sq=False) as f:
        header = f.header.to_dict()
        for r in f.fetch(until_eof=True):
            seen_total += 1
            fn = r.get_tag("fn") if r.has_tag("fn") else None
            if fns is not None and fn not in fns:
                continue
            name = r.query_name
            if name in recs:
                dup += 1
                continue
            tags = {}
            for t, v, ty in r.get_tags(with_value_type=True):
                if hasattr(v, "tobytes"):
                    v = v.tobytes()
                tags[t] = digest((ty, v))
            recs[name] = {
                "flag": r.flag,
                "seq": digest(r.query_sequence),
                "qual": digest(bytes(r.query_qualities) if r.query_qualities is not None else None),
                "tags": tags,
                "parent": r.get_tag("pi") if r.has_tag("pi") else name,
                "rg": r.get_tag("RG") if r.has_tag("RG") else None,
                "fn": fn,
                "has_mm": r.has_tag("MM") and r.has_tag("ML"),
            }
    return header, recs, dup, seen_total


def summary(path, fns=None):
    header, recs, dup, seen = load(path, fns)
    declared_rg = [h["ID"] for h in header.get("RG", [])]
    pg_ids = [h["ID"] for h in header.get("PG", [])]
    used_rg = collections.Counter(v["rg"] for v in recs.values())
    return {
        "path": path,
        "records_scanned": seen,
        "records": len(recs) + dup,
        "unique_reads": len(recs),
        "duplicates": dup,
        "parents": len({v["parent"] for v in recs.values()}),
        "with_mm_ml": sum(v["has_mm"] for v in recs.values()),
        "per_fn": dict(collections.Counter(v["fn"] for v in recs.values())),
        "rg_declared": len(declared_rg),
        "rg_used": len(used_rg),
        "rg_used_undeclared": sorted(str(k) for k in used_rg if k not in set(declared_rg)),
        "pg_ids": pg_ids,
        "pg_ids_unique": len(pg_ids) == len(set(pg_ids)),
        "co_lines": len(header.get("CO", [])),
        "co_first": (header.get("CO") or [""])[0][:200],
    }, recs


def compare(a, b, fns=None):
    sa, ra = summary(a, fns)
    sb, rb = summary(b, fns)
    ka, kb = set(ra), set(rb)
    common = ka & kb
    field_diff = collections.Counter()
    tag_only = collections.Counter()
    n_diff = 0
    examples = []
    for k in sorted(common):
        x, y = ra[k], rb[k]
        d = []
        for fld in ("flag", "seq", "qual"):
            if x[fld] != y[fld]:
                d.append(fld)
        for t in set(x["tags"]) | set(y["tags"]):
            if t not in x["tags"] or t not in y["tags"]:
                tag_only[t] += 1
                d.append("tag:" + t)
            elif x["tags"][t] != y["tags"][t]:
                d.append("tag:" + t)
        if d:
            n_diff += 1
            field_diff.update(d)
            if len(examples) < 5:
                examples.append({"read": k, "fields": d})
    pa = {v["parent"] for v in ra.values()}
    pb = {v["parent"] for v in rb.values()}
    return {
        "a": sa, "b": sb,
        "reads_only_a": len(ka - kb), "reads_only_b": len(kb - ka),
        "reads_only_a_examples": sorted(ka - kb)[:5], "reads_only_b_examples": sorted(kb - ka)[:5],
        "common_reads": len(common),
        "reads_identical": len(common) - n_diff,
        "reads_differing": n_diff,
        "field_differences": dict(field_diff),
        "tag_present_in_one_only": dict(tag_only),
        "difference_examples": examples,
        "parents_only_a": len(pa - pb), "parents_only_b": len(pb - pa), "common_parents": len(pa & pb),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["summary", "compare"])
    ap.add_argument("bams", nargs="+")
    ap.add_argument("--fn", default=None, help="comma-separated pod5 file names (fn tag) to keep")
    a = ap.parse_args()
    fns = set(a.fn.split(",")) if a.fn else None
    if a.mode == "summary":
        out = {b: summary(b, fns)[0] for b in a.bams}
    else:
        if len(a.bams) != 2:
            sys.exit("compare needs two BAMs")
        out = compare(a.bams[0], a.bams[1], fns)
    out["fn_filter"] = sorted(fns) if fns else None
    out["pysam"] = pysam.__version__
    json.dump(out, sys.stdout, indent=1, sort_keys=True)
    print()


if __name__ == "__main__":
    main()
