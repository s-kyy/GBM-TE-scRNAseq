#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500

source ./6.4-extGSEAReports.sh gbm_ge_celltypes_markers_all_avgexp_bygroup/gsea_results2/  gbm_ge >stderr 2>&1
