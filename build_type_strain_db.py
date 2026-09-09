#!/usr/bin/env python3
"""
build_type_strain_db.py — assemble a LOCAL, OFFLINE collection of RefSeq
type-strain proteomes for use with Bakta/Prokka's --proteins flag.

WHY THIS EXISTS
---------------
The cluster has no internet access, but the annotation stage wants to hand
Bakta/Prokka the proteome of the RefSeq *type strain* for whatever species
GTDB-Tk assigned to each genome. So the download has to happen somewhere
else and be copied across.

Run this on a machine WITH internet (your laptop, a login node with egress,
etc.), then rsync the resulting directory to the cluster and point
`type_strain_db_dir:` in pipeline_config.yaml at it.

TWO-PHASE WORKFLOW
------------------
  1. Run the pipeline once. GTDB-Tk classifies everything and the annotation
     stage writes <output>/annotation/missing_species.txt listing every
     species it could not resolve to a local proteome.
  2. Copy that file here, run:
         python3 build_type_strain_db.py \\
             --species-list missing_species.txt \\
             --out-dir ./type_strain_db
     rsync ./type_strain_db to the cluster, set type_strain_db_dir, and
     re-run. Annotation picks up the new proteomes.

  (First time through you can also seed it from a GTDB-Tk summary directly:
       --from-gtdbtk gtdbtk.bac120.summary.tsv --from-gtdbtk gtdbtk.ar53.summary.tsv)

WHAT COUNTS AS A TYPE STRAIN
----------------------------
NCBI's assembly_summary_refseq.txt carries a `relation_to_type_material`
column. Any non-empty value means the assembly derives from type material.
We prefer, in order:
    assembly from type material
    assembly from synonym type material
    assembly from pathotype material
    assembly designated as neotype
    (anything else non-empty)
and within a tier prefer Complete Genome > Chromosome > Scaffold > Contig,
then prefer 'reference genome' / 'representative genome' RefSeq categories.

OUTPUT LAYOUT
-------------
    <out-dir>/
        index.tsv                 species \\t genus \\t accession \\t organism \\t rel_path
        proteomes/<Genus>_<species>.faa

Only stdlib is used, so this runs anywhere with Python 3.8+.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ASSEMBLY_SUMMARY_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
)

# Preference order for relation_to_type_material (lower index == better)
TYPE_PREFERENCE = [
    "assembly from type material",
    "assembly from synonym type material",
    "assembly from pathotype material",
    "assembly designated as neotype",
    "assembly designated as reftype",
]
ASSEMBLY_LEVEL_RANK = {
    "Complete Genome": 0,
    "Chromosome": 1,
    "Scaffold": 2,
    "Contig": 3,
}
REFSEQ_CATEGORY_RANK = {
    "reference genome": 0,
    "representative genome": 1,
    "na": 2,
}


# ---------------------------------------------------------------------------
# Taxon-name normalisation (shared conceptually with _resolve_annotation_plan)
# ---------------------------------------------------------------------------

def normalize_taxon(name: str) -> str:
    """Strip GTDB placeholder suffixes and tidy whitespace.

    GTDB appends capital-letter suffixes to split polyphyletic taxa, e.g.
    'Methanobrevibacter_A smithii' or 'Escherichia coli_D'. NCBI/RefSeq uses
    the unsuffixed name, so we strip them for lookup purposes.
    """
    if not name:
        return ""
    name = name.strip()
    name = re.sub(r"^[a-z]__", "", name)          # drop d__/p__/g__/s__ prefix
    parts = [re.sub(r"_[A-Z]+$", "", p) for p in name.split()]
    return " ".join(p for p in parts if p).strip()


def species_key(species: str) -> str:
    """Canonical lookup key: normalised, lowercase, single-spaced."""
    return " ".join(normalize_taxon(species).lower().split())


def safe_filename(species: str) -> str:
    s = normalize_taxon(species).replace(" ", "_")
    return re.sub(r"[^A-Za-z0-9._-]", "", s)


# ---------------------------------------------------------------------------
# Input collection
# ---------------------------------------------------------------------------

def species_from_gtdbtk(paths: List[Path]) -> List[str]:
    out: List[str] = []
    for p in paths:
        if not p.exists():
            sys.stderr.write(f"warning: {p} not found, skipping\n")
            continue
        with p.open() as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for row in rdr:
                cls = row.get("classification", "") or ""
                for field in cls.split(";"):
                    field = field.strip()
                    if field.startswith("s__"):
                        sp = field[3:].strip()
                        if sp:
                            out.append(sp)
    return out


def species_from_list(path: Path) -> List[str]:
    out: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Tolerate a leading 's__' or a full classification string
            if ";" in line:
                for field in line.split(";"):
                    field = field.strip()
                    if field.startswith("s__"):
                        line = field
                        break
            out.append(line)
    return out


# ---------------------------------------------------------------------------
# assembly_summary handling
# ---------------------------------------------------------------------------

def fetch_assembly_summary(cache: Path, refresh: bool = False) -> Path:
    if cache.exists() and not refresh:
        size_mb = cache.stat().st_size / 1e6
        print(f"[cache] using existing {cache} ({size_mb:.0f} MB)")
        return cache
    cache.parent.mkdir(parents=True, exist_ok=True)
    print(f"[fetch] {ASSEMBLY_SUMMARY_URL}")
    print("        (~500 MB, this takes a few minutes)")
    try:
        with urllib.request.urlopen(ASSEMBLY_SUMMARY_URL, timeout=120) as resp, \
             cache.open("wb") as out:
            copied = 0
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
                copied += len(chunk)
                if copied % (50 << 20) < (1 << 20):
                    print(f"        {copied/1e6:.0f} MB...")
    except urllib.error.URLError as e:
        sys.exit(f"ERROR: could not download assembly summary: {e}")
    print(f"[fetch] saved {cache} ({cache.stat().st_size/1e6:.0f} MB)")
    return cache


def type_rank(rel: str) -> int:
    rel = (rel or "").strip().lower()
    for i, t in enumerate(TYPE_PREFERENCE):
        if rel == t:
            return i
    return len(TYPE_PREFERENCE) if rel else 99


def parse_assembly_summary(path: Path) -> Dict[str, dict]:
    """Return {species_key: best_row} for every species with type material."""
    best: Dict[str, dict] = {}
    n_rows = 0
    n_type = 0
    with path.open(errors="ignore") as fh:
        header: Optional[List[str]] = None
        for line in fh:
            if line.startswith("#"):
                # the last comment line is the header
                if "assembly_accession" in line:
                    header = line.lstrip("#").strip().split("\t")
                continue
            if header is None:
                continue
            n_rows += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue
            row = dict(zip(header, fields))
            rel = row.get("relation_to_type_material", "").strip()
            if not rel:
                continue
            n_type += 1
            organism = row.get("organism_name", "").strip()
            # species = first two tokens of organism_name, typically
            toks = organism.split()
            if len(toks) < 2:
                continue
            sp = f"{toks[0]} {toks[1]}"
            key = species_key(sp)
            if not key:
                continue
            cand = {
                "species": normalize_taxon(sp),
                "genus": normalize_taxon(toks[0]),
                "accession": row.get("assembly_accession", ""),
                "organism": organism,
                "ftp_path": row.get("ftp_path", ""),
                "_type_rank": type_rank(rel),
                "_level_rank": ASSEMBLY_LEVEL_RANK.get(
                    row.get("assembly_level", "").strip(), 9),
                "_cat_rank": REFSEQ_CATEGORY_RANK.get(
                    row.get("refseq_category", "").strip().lower(), 9),
            }
            prev = best.get(key)
            if prev is None or (
                (cand["_type_rank"], cand["_level_rank"], cand["_cat_rank"])
                < (prev["_type_rank"], prev["_level_rank"], prev["_cat_rank"])
            ):
                best[key] = cand
    print(f"[parse] {n_rows:,} assemblies, {n_type:,} from type material, "
          f"{len(best):,} distinct species")
    return best


# ---------------------------------------------------------------------------
# Download proteomes
# ---------------------------------------------------------------------------

def proteome_url(ftp_path: str, kind: str = "faa") -> Optional[str]:
    """URL for the requested artifact.

    kind='gbff' -> <asm>_genomic.gbff.gz   (GenBank; PREFERRED)
    kind='faa'  -> <asm>_protein.faa.gz    (plain protein FASTA)

    GenBank is preferred because Bakta and Prokka both parse it natively and
    pull gene/product/dbxref straight out of the CDS features. Bakta's FASTA
    route instead demands headers shaped '>ID<space>gene~~~product~~~dbxrefs'
    and aborts the whole run on anything else, so GenBank removes an entire
    class of failure.
    """
    if not ftp_path or ftp_path == "na":
        return None
    ftp_path = ftp_path.replace("ftp://", "https://")
    base = ftp_path.rstrip("/").split("/")[-1]
    suffix = "_genomic.gbff.gz" if kind == "gbff" else "_protein.faa.gz"
    return f"{ftp_path}/{base}{suffix}"


def download_proteome(url: str, dest: Path, timeout: int = 120) -> Tuple[bool, str]:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            raw = resp.read()
    except urllib.error.HTTPError as e:
        return False, f"HTTP {e.code}"
    except urllib.error.URLError as e:
        return False, str(e.reason)
    except Exception as e:  # noqa: BLE001 - want to keep going on any failure
        return False, str(e)
    try:
        text = gzip.decompress(raw)
    except OSError as e:
        return False, f"gunzip failed: {e}"
    if not text.strip():
        return False, "empty proteome"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(text)
    return True, ""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="destination directory (rsync this to the cluster)")
    ap.add_argument("--species-list", type=Path, action="append", default=[],
                    help="file with one species per line (repeatable). "
                         "Accepts 's__Genus species' or a full GTDB string.")
    ap.add_argument("--from-gtdbtk", type=Path, action="append", default=[],
                    help="a gtdbtk.*.summary.tsv to harvest species from (repeatable)")
    ap.add_argument("--species", action="append", default=[],
                    help="an explicit species name (repeatable)")
    ap.add_argument("--cache", type=Path, default=None,
                    help="where to keep assembly_summary_refseq.txt "
                         "(default <out-dir>/assembly_summary_refseq.txt)")
    ap.add_argument("--refresh", action="store_true",
                    help="re-download assembly_summary even if cached")
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve species and report, but download nothing")
    ap.add_argument("--format", choices=["gbff", "faa"], default="gbff",
                    help="artifact to fetch for --proteins. 'gbff' (default) "
                         "is GenBank, which Bakta AND Prokka parse natively "
                         "with no header-format pitfalls; 'faa' is plain "
                         "protein FASTA, which the pipeline must reformat "
                         "before Bakta will accept it. Falls back to .faa "
                         "automatically when an assembly has no GenBank file.")
    args = ap.parse_args()

    wanted: List[str] = []
    wanted += args.species
    for p in args.species_list:
        wanted += species_from_list(p)
    if args.from_gtdbtk:
        wanted += species_from_gtdbtk(args.from_gtdbtk)

    # dedupe, keep original spelling for reporting
    seen = set()
    uniq: List[str] = []
    for s in wanted:
        k = species_key(s)
        if k and k not in seen:
            seen.add(k)
            uniq.append(s)

    if not uniq:
        sys.exit("ERROR: no species requested. Use --species-list, "
                 "--from-gtdbtk, or --species.")

    print(f"[input] {len(uniq)} distinct species requested")

    out_dir: Path = args.out_dir
    prot_dir = out_dir / "proteomes"
    out_dir.mkdir(parents=True, exist_ok=True)
    prot_dir.mkdir(parents=True, exist_ok=True)

    cache = args.cache or (out_dir / "assembly_summary_refseq.txt")
    summary = fetch_assembly_summary(cache, refresh=args.refresh)
    type_strains = parse_assembly_summary(summary)

    # Genus-level index for fallback
    by_genus: Dict[str, dict] = {}
    for rec in type_strains.values():
        g = rec["genus"].lower()
        prev = by_genus.get(g)
        if prev is None or (
            (rec["_type_rank"], rec["_level_rank"], rec["_cat_rank"])
            < (prev["_type_rank"], prev["_level_rank"], prev["_cat_rank"])
        ):
            by_genus[g] = rec

    # Merge with any pre-existing index so repeat runs are additive.
    # Entries whose .faa has gone missing (partial rsync, cleaned scratch,
    # interrupted earlier run) are dropped so they get re-fetched rather than
    # silently pointing the pipeline at a file that is not there.
    index_path = out_dir / "index.tsv"
    existing: Dict[str, List[str]] = {}
    stale = 0
    if index_path.exists():
        with index_path.open() as fh:
            rdr = csv.reader(fh, delimiter="\t")
            next(rdr, None)  # header
            for row in rdr:
                if len(row) < 5:
                    continue
                faa = out_dir / row[4]
                if faa.is_file() and faa.stat().st_size > 0:
                    existing[row[0].lower()] = row
                else:
                    stale += 1
        print(f"[index] reusing existing DB at {out_dir}: "
              f"{len(existing)} valid entries"
              + (f", {stale} stale (proteome missing, will re-fetch)" if stale else ""))
    else:
        print(f"[index] no existing DB at {out_dir}, starting fresh")

    resolved, reused, missing, failed = 0, 0, [], []
    for sp in uniq:
        key = species_key(sp)
        if key in existing:
            print(f"  [have] {sp} -> {existing[key][4]}")
            resolved += 1
            reused += 1
            continue
        rec = type_strains.get(key)
        level = "species"
        if rec is None:
            g = normalize_taxon(sp).split()[0].lower() if normalize_taxon(sp) else ""
            rec = by_genus.get(g)
            level = "genus-fallback"
        if rec is None:
            print(f"  [MISS] {sp} — no RefSeq type strain found")
            missing.append(sp)
            continue

        stem = safe_filename(rec["species"])
        kind = "faa" if args.format == "faa" else "gbff"
        ext = ".faa" if kind == "faa" else ".gbff"
        fname = stem + ext
        dest = prot_dir / fname
        rel = f"proteomes/{fname}"

        if dest.exists() and dest.stat().st_size > 0:
            print(f"  [have] {sp} -> {rel} ({level})")
        elif args.dry_run:
            print(f"  [dry ] {sp} -> {rec['accession']} {rec['organism']} "
                  f"({level}, {kind})")
            resolved += 1
            continue
        else:
            url = proteome_url(rec["ftp_path"], kind)
            if not url:
                print(f"  [MISS] {sp} — no ftp_path for {rec['accession']}")
                missing.append(sp)
                continue
            ok, err = download_proteome(url, dest)
            if not ok and kind == "gbff":
                # Not every assembly ships a .gbff; fall back to protein FASTA
                # (the pipeline reformats it for Bakta at annotation time).
                print(f"  [warn] {sp} — no GenBank ({err}), falling back to .faa")
                fname = stem + ".faa"
                dest = prot_dir / fname
                rel = f"proteomes/{fname}"
                url = proteome_url(rec["ftp_path"], "faa")
                ok, err = download_proteome(url, dest) if url else (False, "no ftp_path")
            if not ok:
                print(f"  [FAIL] {sp} — {err}")
                failed.append((sp, err))
                continue
            if dest.suffix == ".faa":
                n_prot = dest.read_text(errors="ignore").count(">")
                print(f"  [ok  ] {sp} -> {rel} ({n_prot:,} proteins, {level})")
            else:
                n_cds = dest.read_text(errors="ignore").count("     CDS  ")
                print(f"  [ok  ] {sp} -> {rel} ({n_cds:,} CDS, GenBank, {level})")

        existing[key] = [
            species_key(sp), rec["genus"], rec["accession"], rec["organism"], rel
        ]
        resolved += 1

    if not args.dry_run:
        with index_path.open("w", newline="") as fh:
            # lineterminator="\n": csv defaults to CRLF, which would leave a
            # carriage return on the last column for any consumer that reads
            # the file with plain line splitting rather than the csv module.
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["species_key", "genus", "accession", "organism", "proteome_rel_path"])
            for key in sorted(existing):
                w.writerow(existing[key])
        print(f"\n[index] wrote {index_path} ({len(existing)} entries)")

    print(f"\nSummary: {resolved} resolved "
          f"({reused} reused from existing DB, {resolved - reused} newly fetched), "
          f"{len(missing)} missing, {len(failed)} failed")
    if missing:
        mp = out_dir / "unresolved_species.txt"
        mp.write_text("\n".join(missing) + "\n")
        print(f"  unresolved species written to {mp}")
        print("  (these have no RefSeq type strain — annotation will simply "
              "run without --proteins for them, which is fine)")
    if failed:
        print("  download failures:")
        for sp, err in failed:
            print(f"    {sp}: {err}")

    print(f"\nNow copy to the cluster and set in pipeline_config.yaml:")
    print(f"  type_strain_db_dir: /path/on/cluster/{out_dir.name}")


if __name__ == "__main__":
    main()
