#!/bin/bash
# Real-data POP-GWAS and labelled-only analysis (revision 2026-10-20; reply to
# Reviewer 1, Comments 1 and 3, Reviewer 3 minor 1): for each of the 24 proteins of
# results/sens_cf_noPRS_noPC/proteins.csv and each mode of step2_realpop.R
# (full, full_jointPERM1, labonly), on the protein's SNP set of the sensitivity analysis.
#   ./submit_realpop.sh              submits the 72-task array
#   ./submit_realpop.sh 3,7,10       resubmits these tasks only
#SBATCH --job-name=realpop
#SBATCH --output=<PROJECT_DIR>/logs/realpop/realpop_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/realpop/realpop_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=56G
#SBATCH --time=0-6:00
#SBATCH --requeue
BASE=<PROJECT_DIR>
mkdir -p $BASE/logs/realpop
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then sbatch --array=${1:-1-72} "$0"; exit 0; fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
t=$((SLURM_ARRAY_TASK_ID - 1)); pi=$((t / 3 + 2)); mi=$((t % 3))
line=$(sed -n "${pi}p" results/sens_cf_noPRS_noPC/proteins.csv)
prot=$(echo "$line" | cut -d, -f2); idx=$(echo "$line" | cut -d, -f3)
case $mi in
  0) export REALPOP_MODE=full ;;
  1) export REALPOP_MODE=full_jointPERM1 ;;
  2) export REALPOP_MODE=labonly ;;
esac
echo "task $SLURM_ARRAY_TASK_ID: $prot (index $idx) $REALPOP_MODE"
[ -s results/realpop_cf_noPRS_noPC/${prot}_${REALPOP_MODE}.rds ] || Rscript code/post_cf3/realpop/step2_realpop.R $idx || exit 1
echo done
