#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G
#SBATCH --array=0-2
cd ./

scriptdir="$(dirname "$0")" # ./gbm_gsc/1_scrna-seq-mapping/2021-06-01_wang_counts_te
cd "$scriptdir"
echo "$scriptdir"

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
./run_cellranger3.sh ../../0_downloads/2021-05-28_wang/SRR10353962 mapRT_SRR10353962 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
./run_cellranger3.sh ../../0_downloads/2021-05-28_wang/SRR10353961 mapRT_SRR10353961 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
./run_cellranger3.sh ../../0_downloads/2021-05-28_wang/SRR10353960 mapRT_SRR10353960 ./ ../../0_downloads/te-ref/refdata_GRCh38-TE 8 15 
fi

