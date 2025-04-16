#!/bin/bash
#SBATCH --time=0:45:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --array=0

module load StdEnv/2020
module load r/4.0.2

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --vanilla ./3.5-filterDeadCells_gbm.R ./gbm_merged_ge_qc_integrated_integrated_umap_clustered_ANNdoublets/merged_ge_qc_integrated_integrated_umap_clustered_ANNdoublets_filtDf_cluster.rds ./gbm_merged_gte_qc_integrated_integrated_umap_clustered_ANNdoublets/merged_gte_qc_integrated_integrated_umap_clustered_ANNdoublets_filtDf_cluster.rds >gbm_ge_gte_filtDC_jan09.out 2>&1
fi
