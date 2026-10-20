#!/bin/bash
# Sensitivity analysis (revision 2026-10-19; reply to Reviewer 1, Comments 4 and 8):
# for each of the 24 proteins of results/sens_cf_noPRS_noPC/proteins.csv and each of
# four settings, step1_sens.R (cross-fitted synthetic phenotype) then step2_sens.R
# (Syn-PALM and Obs-GWAS on the protein's SNP set):
#   _fs20261101, _fs20261102   other fold seeds, 100 predictors
#   _np50, _np200              main folds, 50 or 200 screened continuous predictors
#   ./submit_sens.sh            submits the 96-task array
#SBATCH --job-name=sens
#SBATCH --output=<PROJECT_DIR>/logs/sens/sens_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sens/sens_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=0-5:00
#SBATCH --requeue
BASE=<PROJECT_DIR>
mkdir -p $BASE/logs/sens
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then sbatch --array=1-96 "$0"; exit 0; fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
t=$((SLURM_ARRAY_TASK_ID - 1)); pi=$((t / 4 + 2)); si=$((t % 4))
line=$(sed -n "${pi}p" results/sens_cf_noPRS_noPC/proteins.csv)
prot=$(echo "$line" | cut -d, -f2); idx=$(echo "$line" | cut -d, -f3)
case $si in
  0) export SENS_TAG=_fs20261101 SENS_FOLD_FILE=data/crossfit_folds_K5_seed20261101.rds SENS_NPRED=100 ;;
  1) export SENS_TAG=_fs20261102 SENS_FOLD_FILE=data/crossfit_folds_K5_seed20261102.rds SENS_NPRED=100 ;;
  2) export SENS_TAG=_np50  SENS_FOLD_FILE=data/crossfit_folds_K5.rds SENS_NPRED=50 ;;
  3) export SENS_TAG=_np200 SENS_FOLD_FILE=data/crossfit_folds_K5.rds SENS_NPRED=200 ;;
esac
TAG=_cf_noPRS_noPC$SENS_TAG
D=results/allukb_protein_pQTL_european_discovery/$prot
echo "task $SLURM_ARRAY_TASK_ID: $prot (index $idx) $TAG"
[ -s $D/${prot}_synthetic_phenotype${TAG}.rds ] || Rscript code/post_cf3/sens/step1_sens.R $idx || exit 1
[ -s $D/${prot}${TAG}_scan_subset.rds ] || SYNPALM_S1_TAG=$TAG SENS_SNP_FILE=results/sens_cf_noPRS_noPC/snps_${prot}.txt \
  Rscript code/post_cf3/sens/step2_sens.R $idx 1 || exit 1
echo done
