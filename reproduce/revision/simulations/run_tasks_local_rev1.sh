#!/bin/bash
# Run array tasks of a revision simulation inside one job with P parallel workers
# (for partitions with a per-user job limit, e.g. test; the same seeds and outputs as
# the SLURM array, and finished outputs are skipped by the R scripts, so this can run
# alongside the queued array).
#   sbatch -p test -c 32 --mem=200G -t 0-12:00 run_tasks_local_rev1.sh <kind> <first> <last> [P=$SLURM_CPUS_PER_TASK] [step=1]
#   (also used on hsph with up to ~100 cores and 3 days)
#   kind: tg_t1e (true_geno_typeIerror_calgeno_rev1.R, task t -> replicate (t-1)/40+1, array (t-1)%40+1)
#         tg_power (true_geno_power_rho075_rev1.R, task t -> array t; env RHO, default 0.75)
#         pb_power (power_bias_rho075_rev1.R power, task t -> array t; env GRMR and RHO)
#         pb_bias  (power_bias_rho075_rev1.R bias 0.5, task t -> array t; env RHO)
#         vc09     (variance_components_rev1.R, missing rate 0.9 split as in submit_vc_split_rev1.sh:
#                   task t -> array (t-1)%20+1, replicates 2b+1..2b+2 with b=(t-1)/20;
#                   rho = 0.75 is fitted first (VC_RHO_FIRST; same draws, see the R script))
#   a negative step runs the tasks in reverse order (e.g. <last> <first> with step -1),
#   so that two workers started from opposite ends do not duplicate work
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/local_%x_%j.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/local_%x_%j.err
BASE=<PROJECT_DIR>
KIND=$1; FIRST=$2; LAST=$3; P=${4:-${SLURM_CPUS_PER_TASK:-1}}; STEP=${5:-1}
module load R/4.4.1-fasrc01
MK=<INTEL_ONEAPI_DIR>/mkl/2024.2/lib/intel64
export LD_LIBRARY_PATH=$MK:$LD_LIBRARY_PATH LD_PRELOAD=$MK/libmkl_rt.so
export MKL_INTERFACE_LAYER=GNU,LP64 MKL_THREADING_LAYER=SEQUENTIAL MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd $BASE
run_one() {
  t=$1
  case $KIND in
    tg_t1e)   Rscript simulation_code/rev1/true_geno_typeIerror_calgeno_rev1.R $(( (t - 1) / 40 + 1 )) $(( (t - 1) % 40 + 1 )) 40 > logs/sim_rev1/local_tg_t1e_$t.log 2>&1 ;;
    tg_power) Rscript simulation_code/rev1/true_geno_power_rho075_rev1.R $t "" "" "" ${RHO:-0.75} > logs/sim_rev1/local_tg_power_r${RHO:-0.75}_$t.log 2>&1 ;;
    pb_power) Rscript simulation_code/rev1/power_bias_rho075_rev1.R power $GRMR $t 50 ${RHO:-0.75} > logs/sim_rev1/local_pb_power_${GRMR}_r${RHO:-0.75}_$t.log 2>&1 ;;
    pb_bias)  Rscript simulation_code/rev1/power_bias_rho075_rev1.R bias 0.5 $t 50 ${RHO:-0.75} > logs/sim_rev1/local_pb_bias_r${RHO:-0.75}_$t.log 2>&1 ;;
    vc09)     b=$(( (t - 1) / 20 )); OMP_NUM_THREADS=1 VC_RHO_FIRST=${VC_RHO_FIRST-0.75} VC_RHO_ONLY=${VC_RHO_ONLY-} Rscript simulation_code/rev1/variance_components_rev1.R 0.9 $(( (t - 1) % 20 + 1 )) 10 10000 $(( 2 * b + 1 )) $(( 2 * b + 2 )) > logs/sim_rev1/local_vc09_$t.log 2>&1 ;;
  esac
  echo "task $t exit $?"
}
export -f run_one; export KIND RHO GRMR
seq $FIRST $STEP $LAST | xargs -P $P -I{} bash -c 'run_one {}'
