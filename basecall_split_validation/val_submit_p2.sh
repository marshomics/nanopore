#!/bin/bash
# Phase 2 of rule R5 criterion 5 (run after wave 2 has finished):
#   neg  CPU driver: faults b, c, d on S3neg and the first half of the incremental test on S3inc
#   e    S3neg re-submitted: chunk 01 (truncated in d) is basecalled again, 02/03 skip, merge succeeds
#   f    S3size: the plan records a different size for chunk 01's pod5 file; task 1 must stop before dorado
#   r1   S3 re-submitted after completion: every task and the merge skip
# S3inc (g) and the comparison job follow in val_submit_p2b.sh, after the driver has re-generated
# S3inc (qsub copies a script when it is submitted). Holds are explicit ids (qalter replaces the
# embedded name holds).
set -euo pipefail
V=/ebio/abt3_scratch/jmarsh/software/science/split_build/val
PIN="-q standard.q@node601,standard.q@node602,standard.q@node603,standard.q@node604,standard.q@node605,standard.q@node606,standard.q@node607,standard.q@node608,standard.q@node609"
J=$V/jobs_p2.tsv
[ ! -s $J ] || { echo "$J exists"; exit 4; }
sub() { local name=$1; shift; local id; id=$(qsub -terse -N "$name" "$@"); id=${id%%.*}; printf '%s\t%s\t%s\n' "$name" "$id" "$(date -u +%FT%TZ)" >> $J; echo "$id"; }
cd $V
# The driver copies S3's chunk BAMs, so it waits for the S3 array (AFTER_S3); r1 must not run while
# S3's own merge or the wave 2 comparison runs (AFTER_S3_MERGE), or two merges could write one file.
: "${AFTER_S3:?set AFTER_S3 to the S3 array job id}" "${AFTER_S3_MERGE:?set AFTER_S3_MERGE to the S3 merge and compare job ids}"
NEG=$(sub bsv2_neg -hold_jid $AFTER_S3 neg_driver.sge)
# f: the plan says chunk 01's file has one byte more than it has
awk -F'\t' 'BEGIN { OFS = "\t" } NR > 1 && $1 == "01" { $3 = $3 + 1 } { print }' S3size/submit_scripts/s01_basecall_chunks.tsv > S3size/submit_scripts/plan.tmp
mv S3size/submit_scripts/plan.tmp S3size/submit_scripts/s01_basecall_chunks.tsv
cd $V/S3size/submit_scripts; F=$(sub bsv2_f_size $PIN -t 1-1 s01_basecall.sh)
cd $V/S3neg/submit_scripts
E=$(sub bsv2_e_s01 $PIN -hold_jid $NEG s01_basecall.sh)
EB=$(sub bsv2_e_s01b -hold_jid $E s01b_merge_bam.sh)
cd $V/S3/submit_scripts
R1=$(sub bsv2_r1_s01 $PIN -hold_jid $AFTER_S3_MERGE s01_basecall.sh)
R1B=$(sub bsv2_r1_s01b -hold_jid $R1 s01b_merge_bam.sh)
for pair in "$E $NEG" "$EB $E" "$R1B $R1" "$R1 $AFTER_S3_MERGE" "$NEG $AFTER_S3"; do set -- $pair; qalter -hold_jid $2 $1 > /dev/null; done
cat $J
for id in $(awk -F'\t' '{print $2}' $J); do
  printf '%s %s pred=%s\n' "$id" "$(qstat -j $id | awk '/^job_name/{print $2}')" "$(qstat -j $id | awk '/^jid_predecessor_list:/{print $2}')"
done
