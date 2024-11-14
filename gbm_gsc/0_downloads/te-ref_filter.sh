#!/bin/bash

scriptdir="$(dirname "$0")/te-ref"
cd "$scriptdir"
echo "$scriptdir"

# rename chromosomes from just numbers to having the chr prefix
sed -E 's/^([0-9]+|[XY])\t/chr\1\t/' GRCh38_Ensembl_rmsk_TE.gtf > GRCh38_Ensembl_rmsk_TE_v2.gtf 

# save the contigs names from the gtf file  
awk 'BEGIN{FS="\t"}{print $1}' GRCh38_Ensembl_rmsk_TE_v2.gtf \
  | sort -u > contigs_te

# save the contig names of the ref genome
awk 'BEGIN{FS="\t"}{print $1}'  ../hg38-refdata/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai \ 
  | sort -u > contigs_ge

echo "$(wc -l contigs_*)"
# 194 contigs_ge
# 452 contigs_te
# 646 total

# save inverse matches of col1 of reference genome in col1 of v2 gtf. 
grep -v -f contigs_ge contigs_te > contigs_te_unique
wc -l contigs_te_unique
# 261

# remove unavailable contigs from GTF file. 
grep -v -f contigs_te_unique GRCh38_Ensembl_rmsk_TE_v2.gtf > GRCh38_Ensembl_rmsk_TE_v3.gtf
wc -l *.gtf
#   4693511 GRCh38_Ensembl_rmsk_TE.gtf
#   4693511 GRCh38_Ensembl_rmsk_TE_v2.gtf
#   4532102 GRCh38_Ensembl_rmsk_TE_v3.gtf
#  13919124 total

