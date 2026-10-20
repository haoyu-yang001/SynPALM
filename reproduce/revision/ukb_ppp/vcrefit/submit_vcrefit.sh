#!/bin/bash
# REML refit + subset scan (code/post_cf3/vcrefit/vcrefit.R), one protein per array task;
# task t = row t+1 of results/vcrefit_cf_noPRS_noPC/proteins.csv (column 4 = protein index).
#   bash code/post_cf3/vcrefit/submit_vcrefit.sh [array=1-60] [source=reml]
# source=mom reruns the production MoM through the same code (validation).
# Finished outputs are skipped, so preempted serial_requeue tasks are simply resubmitted.
#SBATCH --job-name=vcrefit
#SBATCH --output=<PROJECT_DIR>/logs/vcrefit/vcrefit_%a_%A.out
#SBATCH --error=<PROJECT_DIR>/logs/vcrefit/vcrefit_%a_%A.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=72G
#SBATCH --time=0-12:00
#SBATCH --requeue
BASE=<PROJECT_DIR>
mkdir -p $BASE/logs/vcrefit
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then sbatch --array=${1:-1-60} --export=ALL,VCREFIT_SOURCE=${2:-reml} "$0"; exit 0; fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
idx=${VCREFIT_INDEX:-$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" results/vcrefit_cf_noPRS_noPC/proteins.csv | cut -d, -f4)}
echo "task $SLURM_ARRAY_TASK_ID: protein index $idx, source ${VCREFIT_SOURCE:-reml}, host $(hostname)"
Rscript code/post_cf3/vcrefit/vcrefit.R $idx || exit 1
echo done
