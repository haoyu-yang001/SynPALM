#!/bin/bash
#SBATCH --job-name=semisyn_truth
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/semisyn/truth.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/semisyn/truth.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=0-04:00
#SBATCH --requeue
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd <PROJECT_DIR>
Rscript simulation_code/rev1/semisyn/marginal_truth_semisyn.R
