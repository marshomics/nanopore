#!/bin/bash
# R5 criterion 2 control: arm S1' (the single job again), generated with S1's config and pinned to the
# same node class. This is the command that was run, as a one-off, at 2026-10-21T07:32Z (job id in
# results/jobs_s1b.tsv); compare_s1b.sge was then submitted with -hold_jid <that id>.
set -e
V=/ebio/abt3_scratch/jmarsh/software/science/split_build/val
SB=/ebio/abt3_scratch/jmarsh/software/science/split_build
PIN=standard.q@node601,standard.q@node602,standard.q@node603,standard.q@node604,standard.q@node605,standard.q@node606,standard.q@node607,standard.q@node608,standard.q@node609
cd $V
[ ! -e S1b ] || { echo "S1b exists"; exit 4; }
/tmp/global2/jmarsh/code/conda_envs/nanopore/bin/python $SB/new_tree/generate_pipeline.py --pod5-dir $V/pod5_test --samplesheet /ebio/abt3_scratch/jmarsh/nanopore/20261012_cm_conjugation/no_sample_id/20261008_1522_P2S-02398-A_PBM28996_eebc9374/samplesheet_template.txt --output-dir $V/S1b --config $V/cfg_S1.yaml > gen_S1b.log 2>&1
cmp -s S1/submit_scripts/s01_basecall.sh <(sed 's#/S1b/#/S1/#g' S1b/submit_scripts/s01_basecall.sh) && echo "S1b s01 = S1 s01 (paths aside)"
cd S1b/submit_scripts && id=$(qsub -terse -N bsv2_S1b_s01 -q $PIN s01_basecall.sh) && echo "S1b=$id pred=$(qstat -j $id | awk '/^jid_predecessor_list:/{print $2}')"
printf 'bsv2_S1b_s01\t%s\t%s\n' $id "$(date -u +%FT%TZ)" >> $V/jobs_s1b.tsv
