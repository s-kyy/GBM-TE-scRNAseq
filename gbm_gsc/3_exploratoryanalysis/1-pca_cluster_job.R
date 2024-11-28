#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=36G
#SBATCH --mail-user=samantha.yuen@umontreal.ca
#SBATCH --array=0-3

module load StdEnv/2020
module load r/4.0.2

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20241122_merged_bhaduriGBM_wangGBM/merged_gte_qc_integrated_umap.rds >gbm_gte_qc_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20241122_merged_bhaduriGBM_wangGBM/merged_ge_qc_integrated.umap.rds >gbm_ge_qc_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/ge_qc_integrated_umap.rds >healthy_ge_qc_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/gte_qc_integrated_umap.rds >healthy_gte_qc_pca.out 2>&1
fi