#!/usr/bin/env python3
"""Score one validation condition: sv_genotype.py output vs the simulated truth.

Read names carry their source genome as a prefix (A_, B_, INV_, DUP_), so the
allele each spanning read was assigned to can be checked read by read.
"""
import argparse
import json
import os
import re

TAG = {"T002_i2": "A", "T002_i2_with_deletion": "B", "T002_i2_inv3k": "INV", "T002_i2_dup5k": "DUP"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, help="engine run with the known genome versions")
    ap.add_argument("--disc-outdir", default=None, help="discovery-only engine run (scored for discovery)")
    ap.add_argument("--cond", required=True, help="condition id")
    ap.add_argument("--pair", required=True, help="DEL, INV or DUP")
    ap.add_argument("--consensus", required=True, help="genome name used as consensus")
    ap.add_argument("--alt", required=True, help="genome name of the alternative version")
    ap.add_argument("--depth", required=True)
    ap.add_argument("--fraction", type=float, required=True, help="true (realised) alt fraction of reads in the mixture")
    ap.add_argument("--nominal-fraction", type=float, default=None)
    ap.add_argument("--event-pos", type=int, required=True, help="1-based event start in consensus coordinates")
    ap.add_argument("--event-len", type=int, required=True)
    a = ap.parse_args()
    summ = json.load(open(os.path.join(a.outdir, "summary.json")))
    loci = summ["samples"]["mix"]["loci"]
    out = {"cond": a.cond, "pair": a.pair, "consensus": a.consensus, "depth": a.depth, "true_alt_fraction": a.fraction,
           "nominal_alt_fraction": a.nominal_fraction}
    # known-variant locus
    known = [x for x in loci if x.get("source") == "known"]
    if len(known) != 1:
        out["known_status"] = f"{len(known)} known loci placed"
    else:
        k = known[0]
        alt = k["alleles"].get(a.alt, {})
        out.update(known_status=k["status"], known_spanning=k["n"], known_alt_reads=alt.get("reads"),
                   known_alt_fraction=alt.get("fraction"), known_ci_lo=(alt.get("ci95") or [None, None])[0],
                   known_ci_hi=(alt.get("ci95") or [None, None])[1], known_floor=k.get("floor"),
                   known_coverage_ratio=k.get("coverage_ratio"), known_window=k.get("window_bp"),
                   known_eligible=k.get("eligible_reads"), known_used=k.get("used_reads"))
        if alt.get("fraction") is not None:
            out["known_abs_error"] = abs(alt["fraction"] - a.fraction)
            out["known_ci_covers_truth"] = alt["ci95"][0] <= a.fraction <= alt["ci95"][1]
        # read-level assignment check
        rdir = os.path.join(a.outdir, "mix", "reads")
        wrong = total = 0
        for allele in k["alleles"]:
            p = os.path.join(rdir, f"{k['locus']}.{allele}.txt")
            if not os.path.exists(p):
                continue
            want = TAG.get(allele)
            for line in open(p):
                nm = line.strip()
                if not nm:
                    continue
                total += 1
                if nm.split("_", 1)[0] != want:
                    wrong += 1
        out.update(known_reads_checked=total, known_reads_misassigned=wrong)
        j = k.get("junction") or {}
        key = {"DEL": "segment_absent", "INV": "inverted", "DUP": "end_to_start_per_genome"}[a.pair]
        est = (j.get("estimates") or {}).get(key)
        out.update(known_j_status=j.get("status"), known_j_call=j.get("call"), known_j_note=j.get("note"),
                   known_j_counts=json.dumps(j.get("counts")) if j.get("counts") else None)
        if est and est[0] is not None:
            out.update(known_j_fraction=est[0], known_j_ci_lo=est[1], known_j_ci_hi=est[2],
                       known_j_abs_error=abs(est[0] - a.fraction), known_j_ci_covers_truth=est[1] <= a.fraction <= est[2])
    merged = summ["samples"]["mix"].get("merged_into", [])
    out["known_run_merged_discovered"] = ";".join(f"{a}>{b}" for a, b in merged)
    # discovery (scored on the discovery-only run when given)
    if a.disc_outdir:
        dsum = json.load(open(os.path.join(a.disc_outdir, "summary.json")))
        loci = dsum["samples"]["mix"]["loci"]
        out["disc_elapsed_s"] = dsum.get("elapsed_s")
        out["disc_floors"] = json.dumps(dsum["samples"]["mix"].get("floors", {}))
    disc = [x for x in loci if x.get("source") == "discovered"]
    hit, fps, notsup = None, [], 0
    for x in disc:
        m = re.match(r"disc_mix_(DEL|INS|INV|DUP)(\d+)_.+_(\d+)$", x["locus"])
        if not m:
            continue
        kind, size, pos = m.group(1), int(m.group(2)), int(m.group(3))
        reported = x["status"] in ("mixture", "consensus_discrepancy")
        span = size if kind in ("DEL", "INV", "DUP") else 0
        near = pos <= a.event_pos + a.event_len + 1000 and pos + span >= a.event_pos - 1000
        if near:
            if hit is None or (reported and hit["status"] not in ("mixture", "consensus_discrepancy")):
                hit = dict(x, _kind=kind, _size=size)
        elif reported:
            fps.append(f"{x['locus']}:{x['status']}")
        else:
            notsup += 1
    if hit:
        alt_frac = (hit["alleles"].get("ALT_total") or hit["alleles"].get("ALT") or {}).get("fraction")
        # the ALT allele of a discovered locus is the version absent from the consensus: for the DEL pair on the
        # deletion-version consensus Sniffles2 reports the segment as an insertion, i.e. ALT = the full genome
        truth = 1 - a.fraction if (a.pair == "DEL" and hit["_kind"] == "INS") else a.fraction
        out.update(disc_found=True, disc_kind=hit["_kind"], disc_size=hit["_size"], disc_status=hit["status"],
                   disc_method=hit.get("method"), disc_alt_fraction=alt_frac, disc_truth=truth)
        if alt_frac is not None:
            out["disc_abs_error"] = abs(alt_frac - truth)
        j = hit.get("junction") or {}
        key = {"DEL": "segment_absent", "INS": "segment_present", "INV": "inverted",
               "DUP": "end_to_start_per_genome"}[hit["_kind"]]
        est = (j.get("estimates") or {}).get(key)
        out.update(disc_j_status=j.get("status"), disc_j_call=j.get("call"))
        if est and est[0] is not None:
            out.update(disc_j_fraction=est[0], disc_j_abs_error=abs(est[0] - truth),
                       disc_j_ci_covers_truth=est[1] <= truth <= est[2])
    else:
        out.update(disc_found=False)
    out.update(disc_candidates=len(disc), disc_false_reported=len(fps), disc_false_list=";".join(fps[:5]),
               disc_not_supported=notsup, floors=json.dumps(summ["samples"]["mix"].get("floors", {})),
               elapsed_s=summ.get("elapsed_s"))
    print(json.dumps(out))


if __name__ == "__main__":
    main()
