#!/bin/bash
#$ -N svmap_INV_on_A
#$ -cwd
#$ -o /ebio/abt3_scratch/jmarsh/software/science/sv_build/maps/logs/INV_on_A.log
#$ -j y
#$ -l h_vmem=3G
#$ -l h_rt=12:0:0
#$ -pe parallel 8
#$ -hold_jid svsim_INV
. ~/.bashrc >/dev/null 2>&1; mamba activate sv
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 TERM=dumb
set -euo pipefail
cd /ebio/abt3_scratch/jmarsh/software/science/sv_build/maps
echo "node $(hostname) slots $NSLOTS start $(date -u +%FT%TZ)"
ls /ebio/abt3_scratch/jmarsh/software/science/sv_build/sim/pools/INV/part_*.fastq.gz > INV.parts.A
echo "parts: $(wc -l < INV.parts.A)"
if [ ! -s INV_on_A.bam ]; then
  minimap2 -ax map-ont -t $NSLOTS /ebio/abt3_scratch/jmarsh/software/science/sv_build/truth/T002_i2.fasta $(cat INV.parts.A) 2> INV_on_A.minimap2.log \
    | samtools sort -@ 4 -m 1G -o INV_on_A.bam.tmp -
  mv INV_on_A.bam.tmp INV_on_A.bam
fi
samtools index INV_on_A.bam
if [ ! -s INV.bases ]; then
  zcat $(cat INV.parts.A) | awk 'NR % 4 == 2 { n += length($0); r++ } END { print n, r }' > INV.bases.tmp && mv INV.bases.tmp INV.bases
fi
echo "bases/reads: $(cat INV.bases)"
samtools flagstat INV_on_A.bam | head -8
echo "end $(date -u +%FT%TZ)"
