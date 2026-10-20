#!/bin/bash
# List replicates 1..N without output; with "resubmit", resubmit them unless
# they are still queued/running:   ./check_semisyn.sh [N=100] [resubmit]
BASE=<PROJECT_DIR>
N=${1:-100}; D=$BASE/simulation_code/results_rev1/semisyn/reps
active=$(squeue -h -u $USER -n semisyn -r -o "%K" 2>/dev/null | sort -n | tr '\n' ' ')
miss=()
for r in $(seq 1 $N); do f=$(printf "%s/rep_%03d.rds" $D $r); [ -s $f ] || miss+=($r); done
echo "done: $((N - ${#miss[@]}))/$N; missing: ${miss[*]}"; echo "active tasks: $active"
if [ "$2" = resubmit ]; then
  todo=(); for r in "${miss[@]}"; do [[ " $active " == *" $r "* ]] || todo+=($r); done
  [ ${#todo[@]} -gt 0 ] && sbatch --array=$(IFS=,; echo "${todo[*]}") --export=ALL,NTREE=300,SEMISYN_SUBDIR=reps $BASE/simulation_code/rev1/semisyn/submit_semisyn.sh
fi
