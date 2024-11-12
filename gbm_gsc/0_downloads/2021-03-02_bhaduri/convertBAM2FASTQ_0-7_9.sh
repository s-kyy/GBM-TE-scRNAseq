#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G
#SBATCH --array=0-8

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11159_1_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_00
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11159_2_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_01
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11209_1_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_02
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11209_2_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_03
fi

if [[ $SLURM_ARRAY_TASK_ID == 4 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11215_1_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_04
fi

if [[ $SLURM_ARRAY_TASK_ID == 5 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11215_2_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_05
fi

if [[ $SLURM_ARRAY_TASK_ID == 6 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11232_1_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_06
fi

if [[ $SLURM_ARRAY_TASK_ID == 7 ]] ; then
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11232_2_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_07
fi

if [[ $SLURM_ARRAY_TASK_ID == 8 ]] ; then 
bamtofastq-1.3.2 --nthreads=8 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/SF11285_2_possorted_genome_bam.bam.1 ~/scratch/GBM-TE-scRNAseq/0_downloads/2021-03-02_bhaduri/data_09 
fi 

