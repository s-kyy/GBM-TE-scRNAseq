#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   
#SBATCH --mem-per-cpu=32G
#SBATCH --array=0-1
#SBATCH --mail-user=samantha.yuen@umontreal.ca
#SBATCH --mail-type=ALL

pwd
module load gcc/9.3.0; module load gdal/3.2.3; module load r/4.0.2

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --no-save --no-restore --verbose /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_ge.R &> r_findallmarkers_ge.Rout
# R CMD BATCH /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_ge.R
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
Rscript --no-save --no-restore --verbose /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_gte.R &> r_findallmarkers_gte.Rout
# R CMD BATCH /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_gte.R
fi

# if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
# # R CMD BATCH /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_ge_0.3.R
# fi

# if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
# R CMD BATCH /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_ge_0.4.R
# fi

# if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
# R CMD BATCH /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_gte_0.3.R
# fi

# if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
# R CMD BATCH /home/samkyy/script/gete-gbm/results/GBMGSCTE/r_findallmarkers_gte_0.4.R
# fi
