#!/bin/sh

# Author: Samantha Yuen
# Created: 2024-08-07
# PURPOSE: Summarize GSEA Results of GBMSC (+TE) PC20r0.3 scRNA-seq
# DESCRIPTION:
# Export GSEA results from .tsv files || Example filepath: results/GBMSC_GSEA-jul09/cluster0_VarheekWang_GE_table.Gsea.1625866268331/gsea_report_for_cluster0_1625866268331.tsv
# INPUT: takes one parameter for GSEA output directory. 
# OUTPUT: .tsv file combining .tsv files for all clusters. Includes mapping information. .tsv files

GSEA_PATH=${1}
OUTPUT_PATH=${2}

echo "creating output folder"
cd $OUTPUT_PATH
mkdir "summary"
cd summary
OUTPUT_PATH="$(pwd)"

#if [ -f temp ] && rm temp; then
#touch temp
#fi

#if [ -f gte_temp ] && rm gte_temp; then
#touch gte_temp
#fi 

echo "moving to gsea path"
cd $GSEA_PATH

# loop through each output folder (store name)

for dir in ${GSEA_PATH}*/     # list directories in the form "/tmp/dirname/"
do
    DIR_PATH=${dir%*/}      # remove the trailing "/"
    DIR_NAME="${DIR_PATH##*/}"    # print everything after the final "/"

    echo "Summarizing folder: ${DIR_PATH}"
    echo "${OUTPUT_PATH}/${DIR_NAME}"

    if [ -f "${OUTPUT_PATH}/${DIRNAME}.tsv" ] && rm "${OUTPUT_PATH}/${DIR_NAME}.tsv"; then 
    	touch "${OUTPUT_PATH}/${DIR_NAME}.tsv"
    fi

    cd ${DIR_PATH}
    IFS=$'\n' clusters=($(ls | grep "_table" | sed 's/.*\(^cluster[[:digit:]]*\).*/\1/'))
    
    # loop through each cluster and output to folder

    for cluster in "${clusters[@]}"; do
        FILE=$(find -path "*_table/gsea_report_for_${i}_*.tsv" -type f)
	echo "${FILE}"
        cut -f1,4- ${FILE} | tail -n +2 | sed "s/^/${DIR_NAME}\t${cluster}\t/" >> "${OUTPUT_PATH}/${DIRNAME}.tsv"

    done # cluster loop end

    HEADER=$(cut -f1,4- ${FILE} | head -n 1)
    sed -i "1s/^/map\tcluster\t${HEADER}\n/" "${OUTPUT_PATH}/${DIRNAME}.tsv"
    
    echo "Extraction of GSEA files ${DIRNAME} complete"

done # dir loop end

echo "See output gseaResults.tsv in ${OUTPUT_PATH}"

# save and end script
