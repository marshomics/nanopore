#!/bin/bash
#$ -N svval
#$ -cwd
#$ -o /ebio/abt3_scratch/jmarsh/software/science/sv_build/val/logs/val.$TASK_ID.log
#$ -j y
#$ -l h_vmem=4G
#$ -l h_rt=8:0:0
#$ -pe parallel 4
#$ -t 1-120
#$ -tc 30
. ~/.bashrc >/dev/null 2>&1; mamba activate sv
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 TERM=dumb
set -euo pipefail
S=/ebio/abt3_scratch/jmarsh/software/science/sv_build
line=$(awk -F'\t' -v t="$SGE_TASK_ID" 'NR == t + 1' "$S/val/conditions.tsv")
IFS=$'\t' read -r cond pair cons_name cons_tag base_pool alt_pool depth frac g_base g_alt alt_name epos elen glen <<< "$line"
OUT="$S/val/results/$cond.json"
[ -s "$OUT" ] && { echo "done already: $cond"; exit 0; }
echo "node $(hostname) slots $NSLOTS cond $cond start $(date -u +%FT%TZ)"
T=$(mktemp -d "$S/val/tmp/$cond.XXXXXX")
cd "$T"
read -r bb _ < "$S/maps/$base_pool.bases"
read -r ab _ < "$S/maps/$alt_pool.bases"
sub() {   # bam fraction seed out ; returns 1 when nothing is taken
  local x=$2
  if awk -v x="$x" 'BEGIN { exit !(x <= 0) }'; then return 1; fi
  if awk -v x="$x" 'BEGIN { exit !(x >= 0.999999) }'; then ln -s "$1" "$4"; return 0; fi
  samtools view -@ 2 -b -s "$3$(awk -v x="$x" 'BEGIN { printf ".%06d", int(x * 1000000 + 0.5) }')" -o "$4" "$1"
}
fb=$(awk -v d="$depth" -v f="$frac" -v b="$bb" -v g="$glen" 'BEGIN { printf "%.6f", d * (1 - f) * g / b }')
fa=$(awk -v d="$depth" -v f="$frac" -v a="$ab" -v g="$glen" 'BEGIN { printf "%.6f", d * f * g / a }')
parts=()
if sub "$S/maps/${base_pool}_on_${cons_tag}.bam" "$fb" 7 b.bam; then parts+=(b.bam); fi
if sub "$S/maps/${alt_pool}_on_${cons_tag}.bam" "$fa" 8 a.bam; then parts+=(a.bam); fi
if [ ${#parts[@]} -eq 2 ]; then samtools merge -f -@ 2 -o mix.bam "${parts[@]}"; else samtools view -b -o mix.bam "${parts[0]}"; fi
samtools index mix.bam
realised=$(samtools view -F 0x904 mix.bam | awk -v b="$base_pool" -v a="$alt_pool" '{ split($1, p, "_"); if (p[1] == a) na++; else if (p[1] == b) nb++ } END { printf "%.6f %d %d", (na + nb ? na / (na + nb) : 0), nb, na }')
echo "subsample base $fb alt $fa -> realised alt fraction / base reads / alt reads: $realised"
printf 'mix\t%s\tmix.bam\n' "$S/truth/$cons_name.fasta" > samples.tsv
python "$S/engine/sv_genotype.py" --samples samples.tsv --outdir out --variant-genomes "$g_base" "$g_alt" --discovery --junction-mode always \
    --threads "$NSLOTS" --null-loci 20 > engine.log 2>&1
python "$S/engine/sv_genotype.py" --samples samples.tsv --outdir out_disc --discovery --junction-mode always \
    --threads "$NSLOTS" --null-loci 20 > engine_disc.log 2>&1
python "$S/engine/sv_eval.py" --outdir out --disc-outdir out_disc --cond "$cond" --pair "$pair" --consensus "$cons_name" --alt "$alt_name" \
    --depth "$depth" --fraction "${realised%% *}" --nominal-fraction "$frac" --event-pos "$epos" --event-len "$elen" > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
cp out/events.tsv "$S/val/events/$cond.events.tsv"; cp out/null_loci.tsv "$S/val/events/$cond.null.tsv"; cp engine.log "$S/val/engine_logs/$cond.log"
cp out_disc/events.tsv "$S/val/events/$cond.disc_events.tsv"; cp out_disc/null_loci.tsv "$S/val/events/$cond.disc_null.tsv"
cp engine_disc.log "$S/val/engine_logs/$cond.disc.log"
cd "$S/val"; rm -rf "$T"
echo "end $(date -u +%FT%TZ)"
