#!/usr/bin/env python3
"""Positive control for compare_bams.py: copy a BAM, dropping its first read and changing one base of the
second read's sequence and its qs tag. compare(copy, original) must report 1 read only in the original
and 1 differing read (seq and tag:qs).  Usage: make_control.py IN.bam OUT.bam"""
import sys
import pysam
src, dst = sys.argv[1], sys.argv[2]
with pysam.AlignmentFile(src, check_sq=False) as f, pysam.AlignmentFile(dst, "wb", template=f) as out:
    for i, r in enumerate(f.fetch(until_eof=True)):
        if i == 0:
            continue
        if i == 1:
            q = r.query_qualities
            s = r.query_sequence
            r.query_sequence = ("A" if s[0] != "A" else "C") + s[1:]
            r.query_qualities = q
            r.set_tag("qs", float(r.get_tag("qs")) + 1.0, value_type="f")
        out.write(r)
