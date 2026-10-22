#!/bin/bash
# Supplementary Data 1-5 of the revision (as listed in the Data availability section of
# manuscript_R1.tex), zipped from the current result files into results/supplementary_data_rev1/.
#   bash code/sh_run/build_supplementary_data.sh
set -e
BASE=<PROJECT_DIR>
T=$BASE/results/table_results_cf_noPRS_noPC
OUT=$BASE/results/supplementary_data_rev1
mkdir -p $OUT; cd $OUT; rm -f Supplementary_Data_*.zip
z() { n=$1; shift; zip -j -q Supplementary_Data_$n.zip "$@"; }
z 1 $T/category_rates.csv $T/disease_group_rates.csv
z 2 $T/trans_hubs/trans_hubs_SynG.csv $T/trans_hubs/trans_hubs_ObsG.csv
z 3 $T/summary_per_protein.csv
z 4 $T/external_validation_loci.csv
z 5 $T/predictors_candidates.csv $T/predictors_selected_by_fold.csv.gz $T/predictors_selection_frequency.csv
ls -la $OUT
