#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=xxx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40

module load CCEnv arch/avx2 StdEnv/2020 r/4.0.2

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

Rscript --vanilla ./8.2-pseudotime.R ./20250117_gbm_gte_filtDf_cluster/gbm_gte_celltypes_cnv_gbm_monocle.rds 0.6 figs_monocle >0321_gbm_gte_monocle_pseudotime.out 2>&1
