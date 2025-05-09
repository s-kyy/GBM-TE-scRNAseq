#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=xxx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --partition=compute

## infercnv uses a significant amount of memory: 
## https://www.biorxiv.org/content/10.1101/2024.12.18.629083v1.full


module load CCEnv arch/avx2 StdEnv/2023  gcc/12.3 gdal/3.9.1 geos/3.12.0
module load jags/4.3.2
module load r/4.4.0

cd ~/scratch/runs/3_exploratoryanalysis
echo "$(pwd)" #gbm_gsc/3_exploratoryanalysis

Rscript --vanilla ./4-cnvanalysis.R ./gbm_ge_filtDf_cluster/gbm_ge_celltypes.rds ../data/refdata-gex-GRCh38-2020-A/genes/gene_annotations.txt int06_gsctypes 6 >4.4_gbm_ge_cnv_gsctypes.out 2>&1
