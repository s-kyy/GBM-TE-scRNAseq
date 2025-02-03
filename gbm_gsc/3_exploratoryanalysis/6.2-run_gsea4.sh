#!/bin/sh

RES=$1      #[expression_dataset] filepath
PHENO=$2    #[phenotype_label] filepath
GENESET=$3  #[geneset] filepath
PERM=$4     #[permute] gene-set or phenotype (GSEAPreranked is always by gene-set)

LABEL=$5    # Report Label

MAX=$6      #[set-max] default = 500
MIN=$7      #[set-min] default = 5
SEED=$8     #[rnd_seed] number >= 149 digits
OUT=$9      #[out] output directory

GSEA4="/home/samkyy/bin/GSEA_4.1.0"
cwd=$(pwd)

cd ${cwd}
echo ${cwd}

${GSEA4}/gsea-cli.sh GSEAPreranked -res ${RES} -cls ${PHENO} -gmx ${GENESET} -collapse No_Collapse -permute ${PERM} -nperm 1000 -scoring_scheme weighted -set_max ${MAX} -set_min ${MIN} -mode Max_probe -norm meandiv -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed ${SEED} -zip_report false -rpt_label ${LABEL} -create_svgs false -out ${OUT}

echo "GSEAPreranked run complete for ${LABEL} in ${OUT}"

# Example: -res ~/data/2021-05-06_gsea/GE_gsea_table.txt -cls ~/data/2021-05-06_gsea/phenotypes.cls#cluster0_vs_REST -gmx ~/data/2021-05-06_gsea/genesets.gmx -collapse No_Collapse -permute gene_set -nperm 1000 -scoring_scheme weighted -set_max 500 -set_min 5 -mode Max_probe -norm meandiv -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -zip_report false -rpt_label cluster0_VarheekWang_GE -create_svgs false -out ~/data/2021-05-06_gsea/results