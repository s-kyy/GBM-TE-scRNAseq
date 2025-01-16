#!/bin/bash
#SBATCH --time=7:00:00
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
Rscript --vanilla ./3.6-reintegrate_recluster.R ./20250103_gbm_merged_ge_qc_integrated_integrated_umap_clustered_ANNdoublets/merged_ge_qc_integrated_integrated_umap_clustered_ANNdoublets_filtDf_cluster_filtDC.rds 0.5 0.6 gbm figs_filtDC_int >0113_gbm_ge_re-integrate.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
Rscript --vanilla ./3.6-reintegrate_recluster.R ./20250103_gbm_merged_gte_qc_integrated_integrated_umap_clustered_ANNdoublets/merged_gte_qc_integrated_integrated_umap_clustered_ANNdoublets_filtDf_cluster_filtDC.rds 0.5 0.6 gbm figs_filtDC_int >0113_gbm_gte_re-integrate.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
Rscript --vanilla ./3.6-reintegrate_recluster.R ./20250103_healthy_ge_qc_integrated_integrated_umap_clustered_ANNdoublets/ge_qc_integrated_integrated_umap_clustered_ANNdoublets_filtDf_cluster.rds 0.5 0.6 healthy figs_filtDC_int >0113_healthy_ge_re-integrate.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
Rscript --vanilla ./3.6-reintegrate_recluster.R ./20250103_healthy_gte_qc_integrated_integrated_umap_clustered_ANNdoublets/gte_qc_integrated_integrated_umap_clustered_ANNdoublets_filtDf_cluster.rds 0.5 0.6 healthy figs_filtDC_int >0113_healthy_gte_re-integrate.out 2>&1
fi
