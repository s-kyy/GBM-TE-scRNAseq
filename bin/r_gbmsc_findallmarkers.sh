#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --mem-per-cpu=12G
#SBATCH --array=0-3
#SBATCH --mail-user=samantha.yuen@umontreal.ca
#SBATCH --mail-type=ALL

pwd

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
R CMD BATCH /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_ge_0.3.R
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
R CMD BATCH /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_ge_0.4.R
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
R CMD BATCH /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_gte_0.3.R
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
R CMD BATCH /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_gte_0.4.R
fi