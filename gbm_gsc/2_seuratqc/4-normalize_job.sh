#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --array=0-3

module load StdEnv/2020
module load r/4.0.2

cd ~/scratch/runs/2_seuratqc
echo "$(pwd)" #gbm_gsc/2_seuratqc

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --vanilla ./4-normalize.R ./20230320_healthy/ge_qc.rds >healthy_ge_qc_normalize.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
Rscript --vanilla ./4-normalize.R ./20230320_healthy/gte_qc.rds >healthy_gte_qc_normalize.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
Rscript --vanilla ./4-normalize.R ./20241122_merged_bhaduriGBM_wangGBM/merged_gte_qc.rds >gbm_gte_qc_normalize.out 2>&1
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
Rscript --vanilla ./4-normalize.R ./20241122_merged_bhaduriGBM_wangGBM/merged_ge_qc.rds >gbm_ge_qc_normalize.out 2>&1
fi