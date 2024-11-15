#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G
#SBATCH --array=0-12
cd ./
if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9262920 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9262922 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9262923 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9262932 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 4 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9262937 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 5 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9262938 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 6 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9262946 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 7 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9262956 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 8 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9264382 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 9 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9264383 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 10 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9264385 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 11 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9264388 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 12 ]] ; then
../run_cellranger6.sh ../../0_downloads/2023-03-06_bhaduri_healthy map_SRR9264389 ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 
fi

