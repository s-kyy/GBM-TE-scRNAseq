#!/bin/bash

echo "$(pwd)"

echo "download Wang2020 glioma stem cell scRNA-seq data"

prefetch SRR10353960
prefetch SRR10353961
prefetch SRR10353962

echo "fasterq-dump conversion to fastq files + rename for CellRanger compatibility"
fasterq-dump --split-files SRR10353960 -O ./SRR10353960 -t ./tmp &> outSRR10353960
mv /SRR10353960/SRR10353960_1.fastq GBM27_S1_L008_R1_001.fastq.gz
mv /SRR10353960/SRR10353960_2.fastq GBM27_S1_L008_R2_001.fastq.gz

fasterq-dump --split-files SRR10353961 -O ./SRR10353961 -t ./tmp &> outSRR10353961
mv /SRR10353961/SRR10353961_1.fastq GBM28_S21_L008_R1_001.fastq.gz
mv /SRR10353961/SRR10353961_2.fastq GBM28_S21_L008_R2_001.fastq.gz

fasterq-dump --split-files SRR10353962 -O ./SRR10353962 -t ./tmp &> outSRR10353962
mv /SRR10353962/SRR10353962_1.fastq GBM29_S22_L008_R1_001.fastq.gz
mv /SRR10353962/SRR10353962_2.fastq GBM29_S22_L008_R2_001.fastq.gz
