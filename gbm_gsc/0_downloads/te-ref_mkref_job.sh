#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G

scriptdir="$(dirname "$0")"
cd "$scriptdir/te-ref"
echo "$scriptdir"

$HOME/bin/cellranger-3.0.2/cellranger mkref --genome=refdata_GRCh38-TE --fasta=../hg38-refdata/refdata-gex-GRCh38-2020-A/fasta/genome.fa --genes=$scriptdir/GRCh38_Ensembl_rmsk_TE_v3.gtf --nthreads=8 --memgb=10
