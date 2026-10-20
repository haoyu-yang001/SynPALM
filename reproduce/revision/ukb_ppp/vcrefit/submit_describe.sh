#!/bin/bash
# Descriptive part of the variance-component check (code/post_cf3/vcrefit/describe_projected.R).
#   bash code/post_cf3/vcrefit/submit_describe.sh
#SBATCH --job-name=vcr_desc
#SBATCH --output=<PROJECT_DIR>/logs/vcrefit/describe_%j.out
#SBATCH --error=<PROJECT_DIR>/logs/vcrefit/describe_%j.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=0-8:00
#SBATCH --requeue
BASE=<PROJECT_DIR>
mkdir -p $BASE/logs/vcrefit
if [ -z "$SLURM_JOB_ID" ]; then sbatch "$0"; exit 0; fi
module load R/4.4.1-fasrc01
cd $BASE
Rscript code/post_cf3/vcrefit/describe_projected.R
