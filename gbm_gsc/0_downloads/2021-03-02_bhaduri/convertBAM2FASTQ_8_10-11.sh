#!/bin/bash 
#SBATCH --time=12:00:00 
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 
#SBATCH --mem-per-cpu=40G
#SBATCH --array=0-2


if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11247_1_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_10
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11247_2_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_11
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11285_1_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_8
