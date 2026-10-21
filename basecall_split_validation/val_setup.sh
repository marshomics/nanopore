#!/bin/bash
# Set up the R5 validation: test pod5 folders (symlinks to the four smallest pod5 files of the run),
# test configs and four generated pipelines. Generation only; nothing is submitted here.
set -euo pipefail
SB=/ebio/abt3_scratch/jmarsh/software/science/split_build
V=$SB/val
NPY=/tmp/global2/jmarsh/code/conda_envs/nanopore/bin/python
GEN=$SB/new_tree/generate_pipeline.py
PROJ=/ebio/abt3_scratch/jmarsh/nanopore/20261012_cm_conjugation/no_sample_id
POD5=$PROJ/20261008_1522_P2S-02398-A_PBM28996_eebc9374/pod5
SHEET=$PROJ/20261008_1522_P2S-02398-A_PBM28996_eebc9374/samplesheet_template.txt
FILES="PBM28996_eebc9374_9cc5eed7_77.pod5 PBM28996_eebc9374_9cc5eed7_80.pod5 PBM28996_eebc9374_9cc5eed7_81.pod5 PBM28996_eebc9374_9cc5eed7_82.pod5"
cd $V
for d in S1 S3 S3neg S3size; do [ ! -e $V/$d ] || { echo "$V/$d exists; refusing to overwrite"; exit 4; }; done
for p in pod5_test pod5_neg pod5_size; do
  rm -rf $p && mkdir $p
  for f in $FILES; do ln -s $POD5/$f $p/$f; done
done
# Test config: the production run config with the shipped split resources (4 slots x 80G, 12G cards,
# 24 h; amendments A1 and A2), so both arms fit the same node class (standard.q node601-609).
sed -e 's/^sge_gpumem: .*/sge_gpumem: "12G"/' \
    -e 's/^sge_h_rt_gpu: .*/sge_h_rt_gpu: "24:0:0"/' \
    -e 's/^sge_pe_parallel: .*/sge_pe_parallel: 4/' \
    -e 's/^sge_h_vmem_gpu: .*/sge_h_vmem_gpu: "80G"/' \
    $SB/gen/cfg_run.yaml > cfg_base.yaml
for k in sge_gpumem sge_h_rt_gpu sge_pe_parallel sge_h_vmem_gpu; do
  [ "$(grep -c "^$k:" cfg_base.yaml)" = 1 ] || { echo "key $k not exactly once"; exit 5; }
done
{ cat cfg_base.yaml; printf '\nbasecall_split_jobs: 1\n'; } > cfg_S1.yaml
split_keys='basecall_split_gpumem: "12G"\nbasecall_split_h_rt: "24:0:0"\nbasecall_split_slots: 4\nbasecall_split_h_vmem: "80G"\n'
{ cat cfg_base.yaml; printf "\nbasecall_split_jobs: 3\n$split_keys"; printf 'basecall_split_keep_chunks: true\n'; } > cfg_S3.yaml
{ cat cfg_base.yaml; printf "\nbasecall_split_jobs: 3\n$split_keys"; printf 'basecall_split_keep_chunks: false\n'; } > cfg_S3neg.yaml
gen() { $NPY $GEN --pod5-dir $V/$2 --samplesheet $SHEET --output-dir $V/$1 --config $V/$3 > $V/gen_$1.log 2>&1 || { tail -5 $V/gen_$1.log; exit 3; }; }
gen S1 pod5_test cfg_S1.yaml
gen S3 pod5_test cfg_S3.yaml
gen S3neg pod5_neg cfg_S3neg.yaml
gen S3size pod5_size cfg_S3neg.yaml
grep -h '^s01_basecall' gen_*.log
cat S3/submit_scripts/s01_basecall_chunks.tsv
for d in S3 S3neg S3size; do cut -f2,3 $d/submit_scripts/s01_basecall_chunks.tsv | md5sum; done
grep -E '^#\$ -(l|pe|t|hold_jid)' S1/submit_scripts/s01_basecall.sh S3/submit_scripts/s01_basecall.sh S3/submit_scripts/s01b_merge_bam.sh
