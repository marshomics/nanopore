#!/bin/bash
#$ -N svsim_DUP
#$ -cwd
#$ -o /ebio/abt3_scratch/jmarsh/software/science/sv_build/sim/logs/sim_DUP.$TASK_ID.log
#$ -j y
#$ -l h_vmem=4G
#$ -l h_rt=12:0:0
#$ -t 1-8
. ~/.bashrc >/dev/null 2>&1; mamba activate sv
set -o pipefail
OUT=/ebio/abt3_scratch/jmarsh/software/science/sv_build/sim/pools/DUP/part_$(printf %03d $SGE_TASK_ID).fastq.gz
mkdir -p "$(dirname "$OUT")"
[ -s "$OUT" ] && { echo "exists: $OUT"; exit 0; }
echo "node $(hostname) task $SGE_TASK_ID start $(date -u +%FT%TZ)"
badread simulate --reference /ebio/abt3_scratch/jmarsh/software/science/sv_build/truth/T002_i2_dup5k.fasta --quantity 50x --seed $((4000 + SGE_TASK_ID)) \
    --length 11000,10000 --identity 98,99.8,2 --error_model nanopore2023 --qscore_model nanopore2023 \
    2> "${OUT%.fastq.gz}.badread.log" \
  | awk -v t=DUP 'NR % 4 == 1 { sub(/^@/, "@" t "_") } { print }' | gzip -1 > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"
echo "rc=$? end $(date -u +%FT%TZ) reads=$(( $(zcat "$OUT" | wc -l) / 4 ))"
