#!/bin/sh

set -euxo pipefail

# Author: Samantha Yuen
# Created: 2024-08-07
# PURPOSE: Summarize GSEA Results of GBMSC (+TE) PC20r0.3 scRNA-seq
# DESCRIPTION:
# Export GSEA results from .tsv files || Example filepath: results/GBMSC_GSEA-jul09/cluster0_VarheekWang_GE_table.Gsea.1625866268331/gsea_report_for_cluster0_1625866268331.tsv
# INPUT: takes one parameter for GSEA output directory. 
# OUTPUT: .tsv file combining .tsv files for all clusters. Includes mapping information. .tsv files
# Example command: source ./extGSEAReports.sh ~/scratch/gsea/job_outputs/ ~/scratch/gsea/

if [ $# -eq 0 ]; then
  echo >&2 "No arguments supplied. Exiting..."
  exit 1
elif [ $# -eq 1 ]; then
  echo >&2 "Only one arguement supplied. Exiting..."
  exit 1
elif [ -d ${1} ]; then
  echo "GSEA Directory exists. Processing ...\n\n"
  GSEA_PATH=$(readlink -f ${1})   # Obtain absolute path
  GSEA_PATH="${GSEA_PATH}/"
else 
  echo >&2 "GSEA Directory does not exist. Exiting..."
  exit 1
fi

# Trim trailing spaces and underscores from experiment name (e.g. gbm_ge)
SAMPLE=$(sed 's/^[_[:space:]]*//; s/[_[:space:]]*$//' <<< "${2}") 

echo "Inputed Path: ${GSEA_PATH}"
echo "Experiment Name: ${SAMPLE}"

# if output tsv file doesn't exist (removing fives an error), then generate tsv file with directory name.
GSEA_DIR_NAME=${GSEA_PATH%*/}                   # remove the trailing "/"
GSEA_DIR_NAME="${GSEA_DIR_NAME##*/}"            # print everything after the final "/"
OUTPUT="${GSEA_PATH}${GSEA_DIR_NAME}.tsv"
if [ -f "${OUTPUT}" ] && rm "${OUTPUT}"; then   # Create and/or overwrite output file of the same name
    echo "Creating Output File: ${OUTPUT}\n\n"
    touch "${OUTPUT}"
fi

# loop through each folder (store name)
for dir in ${GSEA_PATH}*/ ; do      # list directories 

    DIR_PATH=${dir%*/}              # remove the trailing "/"
    DIR_NAME="${DIR_PATH##*/}"      # print everything after the final "/"
    
    # Parse Directory name for celltype, cluster and geneset names
    TYPE=$(echo "${DIR_NAME}" | sed -En "s/^${SAMPLE}_([A-Za-z0-9\.]*)_([A-Za-z0-9\.-]*)_(.*).GseaPreranked(.*)/\1/p")
    CLUSTER=$(echo "${DIR_NAME}" | sed -En "s/^${SAMPLE}_([A-Za-z0-9\.*)_([A-Za-z0-9\.-]*)_(.*).GseaPreranked(.*)/\2/p")
    GENESET=$(echo "${DIR_NAME}" | sed -En "s/^${SAMPLE}_([A-Za-z0-9\.]*)_([A-Za-z0-9\.-]*)_(.*).GseaPreranked(.*)/\3/p")

    echo "Seperated Directory name into: ${TYPE} ${CLUSTER} ${GENESET}\n"
    
    IFS=$'\n' files=($(find "${DIR_PATH}" -type f -name "gsea_report*.tsv"))
    
    # loop through each tsv and output to file 
    for file in "${files[@]}"; do
        if [ -f ${file} ]; then
            # using `head -n 1` here instead of awk with `set -euxo pipeline` will result in faulty exit code 141
            HEADER=$(cut -f1,4- ${file} |  awk 'FNR <= 1') 
        fi
	    echo "Exporting ${file}"
        cut -f1,4- ${file} |    # Skip 2-3 columns in tab-delimited file
            tail -n +2 |        # Skip header
            sed "s/^/${TYPE}\t${CLUSTER}\t${GENESET}\t/" >> "${OUTPUT}"
                                # Insert metadata from experiment and append to output file. 
    done # experiment loop end
done # dir loop end

# Add header
sed -i "1s/^/Type\tCluster\tGenset\t${HEADER}\n/" "${OUTPUT}"

echo "\n\nExtraction of GSEA files ${DIR_NAME} complete"
echo "\n\nSee output gseaResults.tsv in ${GSEA_PATH}"

# undo script exit settings
set +euxo pipefail