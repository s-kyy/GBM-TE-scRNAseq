#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48G
#SBATCH --array=0-3

module load StdEnv/2020
module load r/4.0.2

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --vanilla ./3.2-finddoublets.R ./20241203_gbm_merged_ge_qc_integrated_umap/merged_ge_qc_integrated_umap_clustered.rds gbm >1223_gbm_ge_finddoublets.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
Rscript --vanilla ./3.2-finddoublets.R ./20241203_gbm_merged_gte_qc_integrated_umap/merged_gte_qc_integrated_umap_clustered.rds gbm >1223_gbm_gte_finddoublets.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
Rscript --vanilla ./3.2-finddoublets.R ./20241203_healthy_ge_qc_integrated_umap/ge_qc_integrated_umap_clustered.rds healthy healthy >1223_healthy_ge_finddoublets.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
Rscript --vanilla ./3.2-finddoublets.R ./20241203_healthy_gte_qc_integrated_umap/gte_qc_integrated_umap_clustered.rds healthy healthy >1223_healthy_gte_finddoublets.out 2>&1
fi
