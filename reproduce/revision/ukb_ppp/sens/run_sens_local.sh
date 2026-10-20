#!/bin/bash
# Run tasks of submit_sens.sh inside one allocation, P at a time (4 cores each).
#   sbatch -p test -c 28 --mem=120G -t 0-6:00 code/post_cf3/sens/run_sens_local.sh <first> <last> [P=7]
BASE=<PROJECT_DIR>
cd $BASE
seq $1 $2 | xargs -P ${3:-7} -I{} bash -c 'SLURM_ARRAY_TASK_ID={} SLURM_CPUS_PER_TASK=4 bash code/post_cf3/sens/submit_sens.sh > logs/sens/sens_{}.out 2> logs/sens/sens_{}.err; echo "task {} exit $?"'
