#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --array=0-3

module load StdEnv/2020
module load r/4.0.2

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --vanilla ./3.2-finddoublets.R ./20241224_gbm_merged_ge_qc_integrated_integrated_umap/merged_ge_qc_integrated_integrated_umap_clustered.rds gbm >1224_gbm_ge_finddoublets.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
Rscript --vanilla ./3.2-finddoublets.R ./20241224_gbm_merged_gte_qc_integrated_integrated_umap/merged_gte_qc_integrated_integrated_umap_clustered.rds gbm >1224_gbm_gte_finddoublets.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
Rscript --vanilla ./3.2-finddoublets.R ./20241224_healthy_ge_qc_integrated_integrated_umap/ge_qc_integrated_integrated_umap_clustered.rds healthy >1224_healthy_ge_finddoublets.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
Rscript --vanilla ./3.2-finddoublets.R ./20241224_healthy_gte_qc_integrated_integrated_umap/gte_qc_integrated_integrated_umap_clustered.rds healthy >1224_healthy_gte_finddoublets.out 2>&1
fi
