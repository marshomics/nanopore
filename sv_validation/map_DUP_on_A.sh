#!/bin/bash
#$ -N svmap_DUP_on_A
#$ -cwd
#$ -o /ebio/abt3_scratch/jmarsh/software/science/sv_build/maps/logs/DUP_on_A.log
#$ -j y
#$ -l h_vmem=3G
#$ -l h_rt=12:0:0
#$ -pe parallel 8
#$ -hold_jid svsim_DUP
. ~/.bashrc >/dev/null 2>&1; mamba activate sv
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 TERM=dumb
set -euo pipefail
cd /ebio/abt3_scratch/jmarsh/software/science/sv_build/maps
echo "node $(hostname) slots $NSLOTS start $(date -u +%FT%TZ)"
ls /ebio/abt3_scratch/jmarsh/software/science/sv_build/sim/pools/DUP/part_*.fastq.gz > DUP.parts.A
echo "parts: $(wc -l < DUP.parts.A)"
if [ ! -s DUP_on_A.bam ]; then
  minimap2 -ax map-ont -t $NSLOTS /ebio/abt3_scratch/jmarsh/software/science/sv_build/truth/T002_i2.fasta $(cat DUP.parts.A) 2> DUP_on_A.minimap2.log \
    | samtools sort -@ 4 -m 1G -o DUP_on_A.bam.tmp -
  mv DUP_on_A.bam.tmp DUP_on_A.bam
fi
samtools index DUP_on_A.bam
if [ ! -s DUP.bases ]; then
  zcat $(cat DUP.parts.A) | awk 'NR % 4 == 2 { n += length($0); r++ } END { print n, r }' > DUP.bases.tmp && mv DUP.bases.tmp DUP.bases
fi
echo "bases/reads: $(cat DUP.bases)"
samtools flagstat DUP_on_A.bam | head -8
echo "end $(date -u +%FT%TZ)"
