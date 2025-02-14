#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=48G

## infercnv uses a significant amount of memory: 
## https://www.biorxiv.org/content/10.1101/2024.12.18.629083v1.full


module load StdEnv/2020
module load gcc/9.3.0
module load jags/4.3.2
module load r/4.1.0

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --vanilla ./4-cnvanalysis.R ./20250117_gbm_ge_filtDf_cluster/gbm_ge_celltypes.rds ../data/refdata-gex-GRCh38-2020-A/genes/gene_annotations.txt int06_gsctypes 6 >0205_4.4_gbm_ge_cnv_gsctypes.out 2>&1
fi