#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --account=xxx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --array=0-7

module load CCEnv arch/avx2 StdEnv/2020 
module load r/4.0.2

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --vanilla ./6.0-findmarkers.R ./gbm_ge_filtDf_cluster/gbm_ge_celltypes.rds int06_celltypes celltype_markers >gbm_ge_findmarkers_int06_celltypes.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
Rscript --vanilla ./6.0-findmarkers.R ./gbm_ge_filtDf_cluster/gbm_ge_celltypes.rds int06_gsctypes celltype_markers >gbm_ge_findmarkers_int06_gsctypes.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
Rscript --vanilla ./6.0-findmarkers.R ./gbm_ge_filtDf_cluster/gbm_ge_celltypes.rds gsc celltype_markers >gbm_ge_findmarkers_gsc.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
Rscript --vanilla ./6.0-findmarkers.R ./gbm_ge_filtDf_cluster/gbm_ge_celltypes.rds oRG celltype_markers >gbm_ge_findmarkers_oRG.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 4 ]] ; then
Rscript --vanilla ./6.0-findmarkers.R ./gbm_gte_filtDf_cluster/gbm_gte_celltypes.rds int06_celltypes celltype_markers >gbm_gte_findmarkers_int06_celltypes.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 5 ]] ; then
Rscript --vanilla ./6.0-findmarkers.R ./gbm_gte_filtDf_cluster/gbm_gte_celltypes.rds int06_gsctypes celltype_markers >gbm_gte_findmarkers_int06_gsctypes.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 6 ]] ; then
Rscript --vanilla ./6.0-findmarkers.R ./gbm_gte_filtDf_cluster/gbm_gte_celltypes.rds gsc celltype_markers >gbm_gte_findmarkers_gsc.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 7 ]] ; then
Rscript --vanilla ./6.0-findmarkers.R ./gbm_gte_filtDf_cluster/gbm_gte_celltypes.rds oRG celltype_markers >gbm_gte_findmarkers_oRG.out 2>&1
fi
