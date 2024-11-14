#!/bin/bash

scriptdir="$(dirname "$0")"
cd "$scriptdir"
echo "$scriptdir"

# rename chromosomes from just numbers to having the chr prefix
sed -E 's/^([0-9]+|[XY])\t/chr\1\t/' < GRCh38_Ensembl_rmsk_TE.gtf > GRCh38_Ensembl_rmsk_TE_v2.gtf

# save the contigs names from the gtf file  
awk 'BEGIN{FS="\t"}{print $1}' GRCh38_Ensembl_rmsk_TE_v2.gtf> contigs_te

# save the contig names of the ref genome
#awk 'BEGIN{FS="\t"}{print $1}' ../hg38-refdata/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai > contigs_ge
awk 'BEGIN{FS="\t"}{print $1}' ~/projects/def-ytanaka/samkyy/resources/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai > contigs_ge

# save inverse matches of col1 of reference genome in col1 of v2 gtf. 
#grep -v -f contigs_te contigs_ge > nomatch
#sort -u nomatch > uniq_missing
#head uniq_missing
grep -v -f contigs_te contigs_ge | sort -u > contigs_te_unique

