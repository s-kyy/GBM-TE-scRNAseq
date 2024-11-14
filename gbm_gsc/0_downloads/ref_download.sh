#!/bin/bash

scriptdir="$(dirname "$0")"
cd "$scriptdir"
echo "$scriptdir"

# download human reference genome assembly
mkdir ./hg38-ref
cd ./hg38-ref
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

# download TE GTF annotation file. 
mkdir ./te-ref
cd ./te-ref
wget https://www.dropbox.com/sh/1ppg2e0fbc64bqw/AAC5dKo3gsv9j-6eq5DHFJdha/GRCh38_Ensembl_rmsk_TE.gtf.gz?dl=1
mv GRCh38_Ensembl_rmsk_TE.gtf.gz?dl=1 GRCh38_Ensembl_rmsk_TE.gtf.gz
gzip -d GRCh38_Ensembl_rmsk_TE.gtf.gz 
