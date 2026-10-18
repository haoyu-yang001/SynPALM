#!/bin/bash
# Height ablation, cross-fit with relatedness-cluster folds: one job per missing rate.
#   sbatch submit_height.sh
#SBATCH --job-name=cf3_height
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_height/h_%a.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_height/h_%a.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=2-00:00
#SBATCH --array=0-4
RATES=(0 0.25 0.5 0.75 0.9)
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd <PROJECT_DIR>
Rscript simulation_code/height_ablation_crossfit_clusterfold.R ${RATES[$SLURM_ARRAY_TASK_ID]} 20
