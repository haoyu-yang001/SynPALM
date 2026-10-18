#!/bin/bash
## SLURM driver for code/post_cf3/compare_tau_cf.R (Supplementary figure
## "twoSGWAS", global-null vs SNP-specific variance components).
##   bash code/post_cf3/submit_compare_tau_cf.sh yhatgwas [n_array=20]   # Yhat GWAS, array job
##   bash code/post_cf3/submit_compare_tau_cf.sh top      [n_array=20]   # merge + top 2,000 (short job)
##   bash code/post_cf3/submit_compare_tau_cf.sh refit    [n_array=100] [tasks=1-n_array] [cpus=4]
##        2,000 refits (20 SNPs per task for n_array = 100), SNPs of a task run
##        in parallel on <cpus> forked workers (~7 GB each)
## Finished segments / arrays (and finished SNPs inside an array) are skipped,
## so preempted serial_requeue tasks simply restart. Logs: logs/compare_tau_cf/.
BASE=<PROJECT_DIR>
MODE=$1
LOG=$BASE/logs/compare_tau_cf; mkdir -p $LOG
PART=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
if [ -z "$SLURM_JOB_ID" ]; then
  case $MODE in
    yhatgwas) N=${2:-20};  sbatch -J ct_yhat  -p $PART -c 1 --mem=24G -t 0-6:00 --requeue --array=1-$N \
                -o $LOG/yhat_%a_%A.out -e $LOG/yhat_%a_%A.err "$0" $MODE $N ;;
    top)      N=${2:-20};  sbatch -J ct_top   -p $PART -c 1 --mem=24G -t 0-1:00 --requeue \
                -o $LOG/top_%j.out -e $LOG/top_%j.err "$0" $MODE $N ;;
    refit)    N=${2:-100}; C=${4:-4}; sbatch -J ct_refit -p $PART -c $C --mem=$((8 * C + 8))G -t 0-12:00 --requeue --array=${3:-1-$N} \
                -o $LOG/refit_%a_%A.out -e $LOG/refit_%a_%A.err "$0" $MODE $N ;;
    *) echo "unknown mode $MODE"; exit 1 ;;
  esac
  exit 0
fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
if [ "$MODE" = top ]; then
  Rscript code/post_cf3/compare_tau_cf.R top $2
else
  Rscript code/post_cf3/compare_tau_cf.R $MODE $SLURM_ARRAY_TASK_ID $2
fi
