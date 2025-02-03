#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=xxx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --partition=compute
#SBATCH --array=0

## infercnv uses a significant amount of memory: 
## https://www.biorxiv.org/content/10.1101/2024.12.18.629083v1.full


module load CCEnv arch/avx2 StdEnv/2020 gcc/9.3.0 gdal/3.5.1 geos/3.10.2
module load jags/4.3.2
module load r/4.1.0

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
Rscript --vanilla ./4-cnvanalysis.R ./20250117_gbm_ge_filtDf_cluster/gbm_ge_celltypes.rds ../data/refdata-gex-GRCh38-2020-A/genes/gene_annotations.gtf int06_celltypes 6 >0122_gbm_ge_cnv.out 2>&1
fi