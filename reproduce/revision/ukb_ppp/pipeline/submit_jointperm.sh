#!/bin/bash
# Null calibration (Fig S16): joint permutation of (Y, Yhat) after step 1, for
# the 10 proteins of the submitted Fig S16.   ./submit_jointperm.sh [seed]
#SBATCH --job-name=cf3_jperm
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_jointperm/jp_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_jointperm/jp_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=0-6:00
BASE=<PROJECT_DIR>
PIDX=(393 560 710 932 1381 1592 1615 1764 2049 2715)   # ccl24 chit1 cst3 enpp5 il18r1 lilrb2 lrpap1 msmb pilra tnfrsf10c
mkdir -p $BASE/logs/noPRS_noPC_jointperm
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --export=ALL,SYNPALM_JOINT_PERM_SEED=${1:-1} --array=0-9 "$0"; exit 0
fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
export SYNPALM_S1_TAG=_cf_noPRS_noPC SYNPALM_VALIDATE=0
Rscript code/synpalm_main_code/step2_protein_fast_jointperm.R ${PIDX[$SLURM_ARRAY_TASK_ID]} 10
