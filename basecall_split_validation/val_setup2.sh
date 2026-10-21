#!/bin/bash
# Phase 2 set-up: pipeline S3inc (incremental re-planning test) on its own pod5 folder.
set -euo pipefail
V=/ebio/abt3_scratch/jmarsh/software/science/split_build/val
NPY=/tmp/global2/jmarsh/code/conda_envs/nanopore/bin/python
GEN=/ebio/abt3_scratch/jmarsh/software/science/split_build/new_tree/generate_pipeline.py
POD5=/ebio/abt3_scratch/jmarsh/nanopore/20261012_cm_conjugation/no_sample_id/20261008_1522_P2S-02398-A_PBM28996_eebc9374/pod5
SHEET=/ebio/abt3_scratch/jmarsh/nanopore/20261012_cm_conjugation/no_sample_id/20261008_1522_P2S-02398-A_PBM28996_eebc9374/samplesheet_template.txt
cd $V
[ ! -e S3inc ] || { echo "S3inc exists"; exit 4; }
rm -rf pod5_inc && mkdir pod5_inc
for f in PBM28996_eebc9374_9cc5eed7_77.pod5 PBM28996_eebc9374_9cc5eed7_80.pod5 PBM28996_eebc9374_9cc5eed7_81.pod5 PBM28996_eebc9374_9cc5eed7_82.pod5; do ln -s $POD5/$f pod5_inc/$f; done
$NPY $GEN --pod5-dir $V/pod5_inc --samplesheet $SHEET --output-dir $V/S3inc --config $V/cfg_S3neg.yaml > gen_S3inc.log 2>&1
grep '^s01_basecall' gen_S3inc.log
cmp <(cut -f2,3 S3/submit_scripts/s01_basecall_chunks.tsv) <(cut -f2,3 S3inc/submit_scripts/s01_basecall_chunks.tsv) && echo "S3inc plan = S3 plan"
cmp S3/submit_scripts/s01_basecall_chunks.tsv S3inc/submit_scripts/s01_basecall_chunks.tsv && echo "identical manifests"
