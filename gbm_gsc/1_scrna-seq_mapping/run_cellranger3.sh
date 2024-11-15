#!/bin/sh

INPUT=$1
OUTPUT=$2
DIR=$3
REF=$4
CPU=$5
MEM=$(($6* $CPU *9))
MEM=$((MEM / 10)) #90% of available memory
CELLRANGER="~/bin/cellranger-3.0.2"

cd $DIR
echo $(pwd)
echo "FASTQ Dir: ${INPUT}"

${CELLRANGER}/cellranger count --chemistry=SC3Pv2 --localmem=${MEM} --localcores=${CPU} --id=${OUTPUT} --transcriptome=${REF} --fastqs=${INPUT} --sample=${OUTPUT}