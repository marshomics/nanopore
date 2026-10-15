#!/bin/bash
#$ -N fna_probe
#$ -o /ebio/abt3_scratch/jmarsh/nanopore/no_sample_id/output/stdout/fna_probe.log
#$ -j y
#$ -cwd
#$ -l h_rt=0:10:0
#$ -l h_vmem=2G
#$ -pe parallel 1
#
# Reproduces, ON A COMPUTE NODE, the exact test that s21_mijamp failed:
#     [ -s <output>/bakta/<sample>/<sample>.fna ]
#
# Run it the same way s21 runs (qsub, same queue) and compare the output with
# what you see on the login node. If the login node says OK and this says
# MISSING, the problem is filesystem visibility from the execution host, not
# the pipeline logic.
#
#   qsub diagnose_fna_visibility.sh
#   cat <output>/stdout/fna_probe.log

OUT=/ebio/abt3_scratch/jmarsh/nanopore/no_sample_id/output
SAMPLES="M_wolinii_c1 M_millerae_c1 M_wolinii_c2 M_millerae_c2 B4_no1 B4_no2 B4_no3 B4_no4"

echo "host      : $(hostname)"
echo "user      : $(whoami)"
echo "date      : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "TMPDIR    : ${TMPDIR:-<unset>}"
echo

echo "=== is the output root visible at all? ==="
ls -ld "$OUT" 2>&1
ls -ld "$OUT/bakta" 2>&1
echo

echo "=== mount backing the path ==="
df -h "$OUT" 2>&1 | tail -2
echo

echo "=== the exact test s21 performs, per sample ==="
printf "%-16s %-9s %-12s %-10s %s\n" SAMPLE "[ -s ]" SIZE TYPE PATH
for s in $SAMPLES; do
    f="$OUT/bakta/$s/$s.fna"
    if [ -s "$f" ]; then t="OK"; else t="FAILS"; fi
    sz=$(stat -c%s "$f" 2>/dev/null || echo "-")
    # -L resolves symlinks; %F reports regular file / symbolic link / etc.
    ty=$(stat -Lc%F "$f" 2>/dev/null || echo "unstatable")
    printf "%-16s %-9s %-12s %-10s %s\n" "$s" "$t" "$sz" "$ty" "$f"
done
echo

echo "=== directory listings (are the files really there, from here?) ==="
for s in $SAMPLES; do
    echo "--- $s ---"
    ls -l "$OUT/bakta/$s/" 2>&1 | head -6
done
echo

echo "=== symlink check (a dangling symlink also makes -s false) ==="
for s in $SAMPLES; do
    f="$OUT/bakta/$s/$s.fna"
    if [ -L "$f" ]; then
        echo "$s: SYMLINK -> $(readlink -f "$f" 2>&1)"
    fi
done
echo "(no output above = none are symlinks)"
echo

echo "=== consensus fallback availability ==="
for s in $SAMPLES; do
    f="$OUT/consensus/$s.fasta"
    if [ -s "$f" ]; then echo "  $s: consensus OK ($(stat -c%s "$f") bytes)"
    else echo "  $s: consensus MISSING"; fi
done
