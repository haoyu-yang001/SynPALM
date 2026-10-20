#!/bin/bash
# one-off design of the semi-synthetic simulation:  sbatch submit_prep_semisyn.sh
#SBATCH --job-name=semisyn_prep
#SBATCH --output=<PROJECT_DIR>/logs/sim_rev1/semisyn/prep.out
#SBATCH --error=<PROJECT_DIR>/logs/sim_rev1/semisyn/prep.err
#SBATCH --partition=xlin_ice,xlin,sapphire,hsph,shared,serial_requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=0-02:00
#SBATCH --requeue
module load R/4.4.1-fasrc01
cd <PROJECT_DIR>
Rscript simulation_code/rev1/semisyn/prep_semisyn.R
