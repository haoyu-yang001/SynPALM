#!/bin/bash
# rerun of power / bias simulations (rho = 0.75 first; rho = 0, 0.25, 0.5 added 2026-10-16):
#   ./submit_power_bias_rev1.sh power <GRMr> [first] [last] [n_snps] [rho=0.75]   (200 arrays)
#   ./submit_power_bias_rev1.sh bias 0.5 [first] [last] [n_snps] [rho=0.75]       (1000 arrays)
#SBATCH --job-name=rev1_pb
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/pb_%x_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/pb_%x_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-06:00
#SBATCH --requeue
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --job-name=rev1_${1}_$2_r${6:-0.75} --array=${3:-1}-${4:-200} --export=ALL,MODE=$1,GRMR=$2,NSNP=${5:-50},RHO=${6:-0.75} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
Rscript simulation_code/rev1/power_bias_rho075_rev1.R $MODE $GRMR $SLURM_ARRAY_TASK_ID $NSNP ${RHO:-0.75}
