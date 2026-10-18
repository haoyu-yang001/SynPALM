#!/bin/bash
# Conditional step-2 scan of one protein, with lead cis variant(s) as covariates
# in every null model (lilrb5: the cis-pQTL rs12975366 makes the method-of-moments
# variance components degenerate). The unconditional scan is kept in
# results/cond_cis_old/ and supplies the rows of the conditioning variants.
#   sbatch submit_step2_cond.sh <protein_index> <rsid[,rsid]>
#SBATCH --job-name=cf3_s2cond
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_cond/s2cond_%j.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_cond/s2cond_%j.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=0-8:00
BASE=<PROJECT_DIR>
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
name=$(head -1 data/df_allUKB_eur1.tsv | cut -f$((48 + $1)))
export SYNPALM_S1_TAG=_cf_noPRS_noPC SYNPALM_VALIDATE=0 SYNPALM_COND_SNPS=$2
export SYNPALM_UNCOND_FILE=$BASE/results/cond_cis_old/${name}_cf_noPRS_noPC_scan_allSNPs.rds
Rscript code/synpalm_main_code/step2_protein_fast.R $1 10
