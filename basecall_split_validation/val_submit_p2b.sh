#!/bin/bash
# Phase 2, step 2 (after the CPU driver bsv2_neg has finished and re-generated S3inc):
#   g  S3inc re-submitted as submit_all.sh would: chunks 01-03 skip, new chunk 04 is basecalled, merge
# then the phase 2 comparison job, held on every phase 2 job.
set -euo pipefail
V=/ebio/abt3_scratch/jmarsh/software/science/split_build/val
PIN="-q standard.q@node601,standard.q@node602,standard.q@node603,standard.q@node604,standard.q@node605,standard.q@node606,standard.q@node607,standard.q@node608,standard.q@node609"
J=$V/jobs_p2.tsv
grep -q '^regenerate_rc=0' $V/results/inc_regenerate.log || { echo "S3inc was not re-generated"; exit 5; }
grep -q '^#\$ -t 1-4$' $V/S3inc/submit_scripts/s01_basecall.sh || { echo "S3inc s01 is not a 4-task array"; exit 5; }
grep -q '^N=4$' $V/S3inc/submit_scripts/s01b_merge_bam.sh || { echo "S3inc merge does not expect 4 chunks"; exit 5; }
sub() { local name=$1; shift; local id; id=$(qsub -terse -N "$name" "$@"); id=${id%%.*}; printf '%s\t%s\t%s\n' "$name" "$id" "$(date -u +%FT%TZ)" >> $J; echo "$id"; }
cd $V/S3inc/submit_scripts
G=$(sub bsv2_g_s01 $PIN s01_basecall.sh)
GB=$(sub bsv2_g_s01b -hold_jid $G s01b_merge_bam.sh)
qalter -hold_jid $G $GB > /dev/null
ids=$(awk -F'\t' '$1 != "bsv2_compare2" {print $2}' $J | paste -sd, -)
cd $V
CMP=$(sub bsv2_compare2 -hold_jid $ids compare_p2.sge)
cat $J
for id in $G $GB $CMP; do
  printf '%s %s pred=%s\n' "$id" "$(qstat -j $id | awk '/^job_name/{print $2}')" "$(qstat -j $id | awk '/^jid_predecessor_list:/{print $2}')"
done
