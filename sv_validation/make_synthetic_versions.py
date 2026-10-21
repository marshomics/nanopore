#!/usr/bin/env python3
"""Build the synthetic genome versions used to validate s18b_sv from a base genome.

  <base>_inv3k.fasta  inversion of positions 400,001-403,000 (1-based, inclusive)
  <base>_dup5k.fasta  tandem duplication of 1,400,001-1,405,000; the copy follows position 1,405,000

Usage: make_synthetic_versions.py T002_i2.fasta outdir/
"""
import os
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_one(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                raise SystemExit(f"{path}: more than one sequence")
            name = line[1:].split()[0]
        else:
            seq.append(line.strip())
    return name, "".join(seq)


def write(path, name, seq, comment):
    with open(path, "w") as fh:
        fh.write(f">{name} {comment}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


base_path, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
name, s = read_one(base_path)
stem = os.path.splitext(os.path.basename(base_path))[0]
a, b = 400000, 403000                      # 0-based [a, b)
inv = s[:a] + s[a:b][::-1].translate(COMP) + s[b:]
write(os.path.join(outdir, f"{stem}_inv3k.fasta"), f"{name}_inv3k", inv,
      "circular=true synthetic inversion of 400001-403000")
a, b = 1400000, 1405000
dup = s[:b] + s[a:b] + s[b:]
write(os.path.join(outdir, f"{stem}_dup5k.fasta"), f"{name}_dup5k", dup,
      "circular=true synthetic tandem duplication of 1400001-1405000")
print(len(s), len(inv), len(dup))
