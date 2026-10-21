#!/bin/bash
# Create the Clair3 and analysis environments for the small-variant stages and
# record versions + model checksums. Idempotent: existing envs are reused.
set -uo pipefail
source /ebio/abt3_projects/software/miniconda3/etc/profile.d/conda.sh
ENVS=/tmp/global2/jmarsh/code/conda_envs
DOR=/ebio/abt3_scratch/jmarsh/software/dorado_212/dorado-2.1.2-linux-x64/bin/dorado
log(){ echo "[$(date -u +%FT%TZ)] $*"; }
mkdir -p out
if ! ls $ENVS/clair3/bin/run_clair3* >/dev/null 2>&1; then
  log "creating clair3 env"
  mamba create -y -q -p $ENVS/clair3 -c conda-forge -c bioconda "clair3>=2.0.2" > out/clair3_create.log 2>&1 \
    || { log "clair3 env create FAILED"; tail -30 out/clair3_create.log; }
else log "clair3 env exists"; fi
if [ ! -x $ENVS/snv/bin/vcfdist ]; then
  log "creating snv env"
  mamba create -y -q -p $ENVS/snv -c conda-forge -c bioconda python=3.12 pysam numpy scipy pandas mappy \
    minimap2 samtools bcftools htslib vcfdist > out/snv_create.log 2>&1 \
    || { log "snv env create FAILED"; tail -30 out/snv_create.log; }
else log "snv env exists"; fi
URL=https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/r1041_e82_400bps_sup_v520/
M=$ENVS/clair3_models/r1041_e82_400bps_sup_v520
mkdir -p $M
files=$(curl -fsSL "$URL" | grep -o 'href="[^"?/]*"' | sed 's/href="//; s/"$//' | grep -v '^\.\.' | sort -u)
log "model index: $(echo $files | tr '\n' ' ')"
for f in $files; do [ -s "$M/$f" ] || curl -fsSL -o "$M/$f" "$URL$f" || log "download $f FAILED"; done
ls -la $M
# ---- record versions (evidence for gate G1)
CL=$(ls $ENVS/clair3/bin/run_clair3.sh 2>/dev/null || ls $ENVS/clair3/bin/run_clair3.py 2>/dev/null | head -1)
{
  echo "dorado=$($DOR --version 2>&1 | tail -1)"
  echo "clair3_entry=$CL"
  echo "clair3=$($ENVS/clair3/bin/python $(dirname $CL)/run_clair3.py --version 2>&1 | tail -1 || true)"
  echo "clair3_sh=$($CL --version 2>&1 | tail -1 || true)"
  echo "clair3_conda=$(conda list -p $ENVS/clair3 2>/dev/null | awk '$1=="clair3"{print $2}')"
  echo "torch=$(conda list -p $ENVS/clair3 2>/dev/null | awk '$1=="pytorch"||$1=="torch"{print $2}' | head -1)"
  echo "vcfdist=$($ENVS/snv/bin/vcfdist --version 2>&1 | head -1)"
  echo "bcftools=$($ENVS/snv/bin/bcftools --version 2>&1 | head -1)"
  echo "minimap2=$($ENVS/snv/bin/minimap2 --version 2>&1 | head -1)"
  echo "samtools=$($ENVS/snv/bin/samtools --version 2>&1 | head -1)"
  echo "python_pkgs=$($ENVS/snv/bin/python -c 'import pysam,scipy,numpy,pandas,mappy; print("pysam",pysam.__version__,"scipy",scipy.__version__,"numpy",numpy.__version__,"pandas",pandas.__version__,"mappy",getattr(mappy,"__version__","?"))' 2>&1 | tail -1)"
  for f in $M/*; do echo "model_file=$(basename $f) $(stat -c %s $f) $(sha256sum $f | cut -d' ' -f1)"; done
} > out/versions.txt
cat out/versions.txt
