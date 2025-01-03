#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48G
#SBATCH --array=0-3

module load StdEnv/2020
module load r/4.0.2

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ./20241231_gbm_merged_ge_qc_integrated_integrated_umap_clustered_doublets/merged_ge_qc_integrated_integrated_umap_clustered_ANNdoublets.rds ./20241231_gbm_merged_gte_qc_integrated_integrated_umap_clustered_doublets/merged_gte_qc_integrated_integrated_umap_clustered_ANNdoublets.rds gbm >0102_gbm_ge_filtDf.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ./20241231_gbm_merged_gte_qc_integrated_integrated_umap_clustered_doublets/merged_gte_qc_integrated_integrated_umap_clustered_ANNdoublets.rds ./20241231_gbm_merged_ge_qc_integrated_integrated_umap_clustered_doublets/merged_ge_qc_integrated_integrated_umap_clustered_ANNdoublets.rds gbm >0102_gbm_gte_filtDf.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ./20241231_healthy_ge_qc_integrated_integrated_umap_clustered_doublets/ge_qc_integrated_integrated_umap_clustered_ANNdoublets.rds healthy ./20241231_healthy_gte_qc_integrated_integrated_umap_clustered_doublets/gte_qc_integrated_integrated_umap_clustered_ANNdoublets.rds >0102_healthy_ge_filtDf.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ./20241231_healthy_gte_qc_integrated_integrated_umap_clustered_doublets/gte_qc_integrated_integrated_umap_clustered_ANNdoublets.rds ./20241231_healthy_ge_qc_integrated_integrated_umap_clustered_doublets/ge_qc_integrated_integrated_umap_clustered_ANNdoublets.rds healthy healthy >0102_healthy_gte_filtDf.out 2>&1
fi
