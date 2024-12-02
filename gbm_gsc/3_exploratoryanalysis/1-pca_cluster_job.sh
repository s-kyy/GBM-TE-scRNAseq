#!/bin/bash
#SBATCH --time=7:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=36G
#SBATCH --array=0-7

module load StdEnv/2020
module load r/4.0.2

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20241122_merged_bhaduriGBM_wangGBM/merged_ge_qc_integrated_regressCC.rds >gbm_ge_qc_rcc_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20241122_merged_bhaduriGBM_wangGBM/merged_gte_qc_integrated_regressCC.rds >gbm_gte_qc_rcc_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20241122_merged_bhaduriGBM_wangGBM/merged_ge_qc_integrated_regressCCdiff.rds >gbm_ge_qc_rccdiff_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20241122_merged_bhaduriGBM_wangGBM/merged_gte_qc_integrated_regressCCdiff.rds >gbm_gte_qc_rccdiff_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 4 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/ge_qc_integrated_regressCC.rds >healthy_ge_qc_rcc_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 5 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/gte_qc_integrated_regressCC.rds >healthy_gte_qc_rcc_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 6 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/ge_qc_integrated_regressCCdiff.rds >healthy_ge_qc_rccdiff_pca.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 7 ]] ; then
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/gte_qc_integrated_regressCCdiff.rds >healthy_gte_qc_rccdiff_pca.out 2>&1
fi
