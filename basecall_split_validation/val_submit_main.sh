#!/bin/bash
# Submit arms S1 and S3 (GPU jobs pinned to standard.q node601-609), the S3 merge, demux of both
# arms and arm P. Job names start with bsv_ and every hold is an explicit job id, so nothing here can
# attach to, or wait for, the jobs of another pipeline run (stages hold on job names).
set -euo pipefail
V=/ebio/abt3_scratch/jmarsh/software/science/split_build/val
PIN="-q standard.q@node601,standard.q@node602,standard.q@node603,standard.q@node604,standard.q@node605,standard.q@node606,standard.q@node607,standard.q@node608,standard.q@node609"
mkdir -p $V/logs
[ ! -s $V/jobs.tsv ] || { echo "jobs.tsv exists; refusing to resubmit"; exit 4; }
P=${PREFIX:-bsv}
sub() { local name=$1; shift; local id; id=$(qsub -terse -N "$name" "$@"); id=${id%%.*}; printf '%s\t%s\t%s\n' "$name" "$id" "$(date -u +%FT%TZ)" >> $V/jobs.tsv; echo "$id"; }
cd $V/S1/submit_scripts
S1=$(sub ${P}_S1_s01 $PIN s01_basecall.sh)
cd $V/S3/submit_scripts
S3=$(sub ${P}_S3_s01 $PIN s01_basecall.sh)
S3B=$(sub ${P}_S3_s01b -hold_jid $S3 s01b_merge_bam.sh)
cd $V/S1/submit_scripts
sub ${P}_S1_s02 -hold_jid $S1 s02_demux.sh > /dev/null
sub ${P}_S1_s02b -hold_jid $S1 s02b_demux_bam.sh > /dev/null
cd $V/S3/submit_scripts
sub ${P}_S3_s02 -hold_jid $S3B s02_demux.sh > /dev/null
sub ${P}_S3_s02b -hold_jid $S3B s02b_demux_bam.sh > /dev/null
cd $V
[ -s $V/P/P.bam ] || sub ${P}_prod prod_extract.sge > /dev/null
# qsub -hold_jid ADDS to a script's embedded name hold; replace the lists with the ids alone
qalter -hold_jid $S3 $S3B > /dev/null
for n in ${P}_S1_s02 ${P}_S1_s02b; do qalter -hold_jid $S1 $(awk -F'\t' -v n=$n '$1 == n {print $2}' $V/jobs.tsv) > /dev/null; done
for n in ${P}_S3_s02 ${P}_S3_s02b; do qalter -hold_jid $S3B $(awk -F'\t' -v n=$n '$1 == n {print $2}' $V/jobs.tsv) > /dev/null; done
cat $V/jobs.tsv
# every hold must be exactly the id given above
for id in $(awk -F'\t' '{print $2}' $V/jobs.tsv); do
  printf '%s %s pred=%s\n' "$id" "$(qstat -j $id | awk '/^job_name/{print $2}')" "$(qstat -j $id | awk '/^jid_predecessor_list:/{print $2}')"
done
