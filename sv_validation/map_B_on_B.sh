#!/bin/bash
#$ -N svmap_B_on_B
#$ -cwd
#$ -o /ebio/abt3_scratch/jmarsh/software/science/sv_build/maps/logs/B_on_B.log
#$ -j y
#$ -l h_vmem=3G
#$ -l h_rt=12:0:0
#$ -pe parallel 8
#$ -hold_jid svsim_B,svmap_B_on_A
. ~/.bashrc >/dev/null 2>&1; mamba activate sv
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 TERM=dumb
set -euo pipefail
cd /ebio/abt3_scratch/jmarsh/software/science/sv_build/maps
echo "node $(hostname) slots $NSLOTS start $(date -u +%FT%TZ)"
ls /ebio/abt3_scratch/jmarsh/software/science/sv_build/sim/pools/B/part_*.fastq.gz > B.parts.B
echo "parts: $(wc -l < B.parts.B)"
if [ ! -s B_on_B.bam ]; then
  minimap2 -ax map-ont -t $NSLOTS /ebio/abt3_scratch/jmarsh/software/science/sv_build/truth/T002_i2_with_deletion.fasta $(cat B.parts.B) 2> B_on_B.minimap2.log \
    | samtools sort -@ 4 -m 1G -o B_on_B.bam.tmp -
  mv B_on_B.bam.tmp B_on_B.bam
fi
samtools index B_on_B.bam
if [ ! -s B.bases ]; then
  zcat $(cat B.parts.B) | awk 'NR % 4 == 2 { n += length($0); r++ } END { print n, r }' > B.bases.tmp && mv B.bases.tmp B.bases
fi
echo "bases/reads: $(cat B.bases)"
samtools flagstat B_on_B.bam | head -8
echo "end $(date -u +%FT%TZ)"
