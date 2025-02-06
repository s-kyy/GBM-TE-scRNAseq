#!/bin/bash

WDIR="$(pwd)" # refdata-gex-GRCh38-2020-A/genes/
GTF=$1 # refdata-gex-GRCh38-2020-A/genes/genes.gtf

# Columns: Gene-name chromosome start-base end-base (tab-delimited)
# Col-IDs: 

echo "${WDIR}"

# export features of given gene to file 
awk -F'[\t; "]' '$3 == "gene"{printf "%s\t%s\t%s\t%s\n", $26, $1, $4, $5}' ${GTF} | awk '!seen[$1]++' > gene_annotations.txt

# Number of genes before removing duplicate lines
# awk -F'[\t; "]' '$3 == "gene"{printf "%s\t%s\t%s\t%s\n", $26, $1, $4, $5}' ${GTF} > gene_annotations.gtf
# wc -l gene_annotations.gtf 
# 36601 gene_annotations.gtf

# Number of genes before after removing duplicate lines
# awk -F'[\t; "]' '$3 == "gene"{printf "%s\t%s\t%s\t%s\n", $26, $1, $4, $5}' ${GTF} > gene_annotations.gtf
# awk '!seen[$1]++' gene_annotations.gtf | wc -l
# 36591