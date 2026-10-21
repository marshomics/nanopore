#!/bin/bash
#$ -N svmap_A_on_B
#$ -cwd
#$ -o /ebio/abt3_scratch/jmarsh/software/science/sv_build/maps/logs/A_on_B.log
#$ -j y
#$ -l h_vmem=3G
#$ -l h_rt=12:0:0
#$ -pe parallel 8
#$ -hold_jid svsim_A,svmap_A_on_A
. ~/.bashrc >/dev/null 2>&1; mamba activate sv
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 TERM=dumb
set -euo pipefail
cd /ebio/abt3_scratch/jmarsh/software/science/sv_build/maps
echo "node $(hostname) slots $NSLOTS start $(date -u +%FT%TZ)"
ls /ebio/abt3_scratch/jmarsh/software/science/sv_build/sim/pools/A/part_*.fastq.gz > A.parts.B
echo "parts: $(wc -l < A.parts.B)"
if [ ! -s A_on_B.bam ]; then
  minimap2 -ax map-ont -t $NSLOTS /ebio/abt3_scratch/jmarsh/software/science/sv_build/truth/T002_i2_with_deletion.fasta $(cat A.parts.B) 2> A_on_B.minimap2.log \
    | samtools sort -@ 4 -m 1G -o A_on_B.bam.tmp -
  mv A_on_B.bam.tmp A_on_B.bam
fi
samtools index A_on_B.bam
if [ ! -s A.bases ]; then
  zcat $(cat A.parts.B) | awk 'NR % 4 == 2 { n += length($0); r++ } END { print n, r }' > A.bases.tmp && mv A.bases.tmp A.bases
fi
echo "bases/reads: $(cat A.bases)"
samtools flagstat A_on_B.bam | head -8
echo "end $(date -u +%FT%TZ)"
