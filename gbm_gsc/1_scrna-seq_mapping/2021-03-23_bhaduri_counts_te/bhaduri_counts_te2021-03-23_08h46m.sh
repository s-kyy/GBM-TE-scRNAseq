#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G
#SBATCH --array=0-11

scriptdir="$(dirname "$0")" #./gbm_gsc/1_scrna-seq-mapping/2021-03-23_bhaduri_counts_te
cd "$scriptdir"
echo "$scriptdir"

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_06/SF11232_1_MissingLibrary_1_HKMVGBCXY mapRT_0 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_05/SF22215_GBM2_MissingLibrary_1_CAUFRANXX mapRT_1 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_01/SF11159_2_MissingLibrary_1_HFMLTBCXY mapRT_2 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_8/Tumor2_1_MissingLibrary_1_CB64GANXX mapRT_3 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 4 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_07/SF11232_2_MissingLibrary_1_HKMVGBCXY mapRT_4 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 5 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_02/SF11209_GBM1_MissingLibrary_1_CAUFRANXX mapRT_5 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 6 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_03/SF11209_GBM2_MissingLibrary_1_CAUFRANXX mapRT_6 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 7 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_00/SF11159_1_MissingLibrary_1_HFMLTBCXY mapRT_7 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 8 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_04/SF22215_GBM1_MissingLibrary_1_CAUFRANXX mapRT_8 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 9 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_10/TQ1_MissingLibrary_1_CAT6PANXX mapRT_9 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 10 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_11/TQ2_MissingLibrary_1_CAT6PANXX mapRT_10 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 11 ]] ; then
../run_cellranger3.sh ../../0_downloads/2021-03-23_bhaduri/data_09/Tumor2_2_MissingLibrary_1_CB64GANXX mapRT_11 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

