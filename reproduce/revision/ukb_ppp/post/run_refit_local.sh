#!/bin/bash
## Run the LEP SNP-specific refits (compare_tau_cf.R refit <a> 100) for arrays
## <first>..<last> inside one job with P single-core workers (for the test partition);
## finished arrays / SNPs are skipped, so it can run alongside the ct_refit array.
##   sbatch -p test -c 24 --mem=215G -t 0-12:00 code/post_cf3/run_refit_local.sh 3 100 24
#SBATCH -J ct_refit_T
#SBATCH -o <PROJECT_DIR>/logs/compare_tau_cf/local_%j.out
#SBATCH -e <PROJECT_DIR>/logs/compare_tau_cf/local_%j.err
BASE=<PROJECT_DIR>
FIRST=$1; LAST=$2; P=${3:-${SLURM_CPUS_PER_TASK:-1}}
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
run_one() { SLURM_CPUS_PER_TASK=1 Rscript code/post_cf3/compare_tau_cf.R refit $1 100 > logs/compare_tau_cf/refit_local_$1.out 2>&1; echo "array $1 exit $?"; }
export -f run_one
seq $FIRST $LAST | xargs -P $P -I{} bash -c 'run_one {}'
