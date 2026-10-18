#!/bin/bash
# SNPsea (Slowikowski et al. 2014; v1.0.4 linux binary, data SNPsea_data_20140520) with
# the Gene Atlas 2004 expression matrix, up to 1e6 null SNP sets per tissue (P resolved
# down to 1e-6, far below the Bonferroni threshold 0.05/79) and 100 exceedances to stop, for the SNP sets of
# code/post_cf3/snpsea_prep.R. Tool and data live outside the project in ~/tools/.
#   sbatch -p test -c 32 --mem=64G -t 0-8:00 run_snpsea.sh [sets, default "SynG ObsG SynOnly"]
# (the threads are shared equally among the sets run in the job)
#SBATCH --job-name=snpsea
#SBATCH --output=<PROJECT_DIR>/logs/noPRS_noPC_post/snpsea_%j.out
#SBATCH --error=<PROJECT_DIR>/logs/noPRS_noPC_post/snpsea_%j.err
BASE=<PROJECT_DIR>
SNPSEA=$HOME/tools/snpsea/bin/snpsea-linux64
DATA=$HOME/tools/snpsea_data
OUT=$BASE/results/table_results_cf_noPRS_noPC/snpsea
OUTRUN=${SNPSEA_OUT:-$OUT}   # output folder of this run (inputs always from $OUT)
mkdir -p $OUTRUN
SETS=${1:-SynG ObsG SynOnly}
T=$(( ${SLURM_CPUS_PER_TASK:-3} / $(echo $SETS | wc -w) ))
for s in $SETS; do
  [ "$(wc -l < $OUTRUN/$s/condition_pvalues.txt 2>/dev/null || echo 0)" -ge 80 ] && continue   # all 79 tissues done
  $SNPSEA --snps $OUT/$s.gwas --gene-matrix $DATA/GeneAtlas2004.gct.gz --gene-intervals $DATA/NCBIgenes2013.bed.gz \
    --snp-intervals $DATA/TGP2011.bed.gz --null-snps $DATA/Lango2010.txt.gz --out $OUTRUN/$s \
    --slop 10e3 --threads $T --null-snpsets 0 --min-observations 100 --max-iterations 1e6 > $OUTRUN/$s.log 2>&1 &
done
wait
ls -la $OUTRUN/*/condition_pvalues.txt
