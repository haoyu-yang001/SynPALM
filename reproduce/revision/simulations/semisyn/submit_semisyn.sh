#!/bin/bash
# Semi-synthetic simulation, one replicate per array task (outputs are skipped
# if present, so requeued / resubmitted tasks are safe):
#   ./submit_semisyn.sh <first> <last> [ntree=300] [subdir=reps]
#SBATCH --job-name=semisyn
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/semisyn/rep_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/semisyn/rep_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=0-08:00
#SBATCH --requeue
BASE=<PROJECT_DIR>
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  sbatch --array=${1:-1}-${2:-1} --export=ALL,NTREE=${3:-300},SEMISYN_SUBDIR=${4:-reps} "$0"; exit 0
fi
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
Rscript simulation_code/rev1/semisyn/run_semisyn_rep.R $SLURM_ARRAY_TASK_ID ${NTREE:-300}
