#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --mem-per-cpu=8G
#SBATCH --array=0-1
#SBATCH --mail-user=samantha.yuen@umontreal.ca
#SBATCH --mail-type=ALL

pwd

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
R CMD BATCH /home/samkyy/scratch/gete-gbm/results/GBMGSCTE/r_findallmarkers_ge.R
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
R CMD BATCH /home/samkyy/scratch/gete-gbm/results/GBMGSCTE/r_findallmarkers_gte.R
fi