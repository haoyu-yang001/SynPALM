#!/bin/bash
#SBATCH --job-name=cf3_vcdiag
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_vcfix/d_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_vcfix/d_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=0-3:00
#SBATCH --array=0-4
PIDX=(1615 560 1381 1764 710)
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd <PROJECT_DIR>
export SYNPALM_S1_TAG=_cf_noPRS_noPC SYNPALM_JOINT_PERM_SEED=1
Rscript code/fastscan/diag_vcfix_lambda.R ${PIDX[$SLURM_ARRAY_TASK_ID]}
