#!/bin/bash
# Variance components MoM / REML / BOLT-REML: ./submit_vc_rev1.sh <miss> [first] [last] [n_rep] [M]
#SBATCH --job-name=rev1_vc
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/vc_%x_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/vc_%x_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=3-00:00
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --job-name=rev1_vc_$1 --array=${2:-1}-${3:-20} --export=ALL,MISS=$1,NREP=${4:-10},MSNP=${5:-10000} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
export OMP_NUM_THREADS=1
cd $BASE
Rscript simulation_code/rev1/variance_components_rev1.R $MISS $SLURM_ARRAY_TASK_ID $NREP $MSNP
