#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=15G

Rscript ~/scratch/gete-gbm/bin/meanofratios.R "~/scratch/gete-gbm/results/2021-07-08/gte_gbm.sc.int.umap0.7_mgmt.RData" &> out
