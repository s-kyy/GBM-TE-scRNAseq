#!/bin/bash
#SBATCH --time=24:00:00
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
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20241122_merged_bhaduriGBM_wangGBM/20241203_gbm_int_umap/20241224_gbm_ge_dimred/merged_ge_qc_integrated_integrated_umap.rds gbm >1227_gbm_ge_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20241122_merged_bhaduriGBM_wangGBM/20241203_gbm_int_umap/20241224_gbm_gte_dimred/merged_gte_qc_integrated_integrated_umap.rds gbm >1227_gbm_gte_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/20241203_healthy_int_umap/20241224_healthy_ge_dimred/ge_qc_integrated_integrated_umap.rds healthy >1227_healthy_ge_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/20241203_healthy_int_umap/20241224_healthy_gte_dimred/gte_qc_integrated_integrated_umap.rds healthy >1227_healthy_gte_pca.out 2>&1
fi
