#!/bin/sh

RES=$1      #[expression_dataset] filepath
PHENO=$2    #[phenotype_label] filepath
GENESET=$3  #[geneset] filepath
PERM=$4     #[permute] gene-set or phenotype

LABEL=$5    # Report Label

MAX=$6      #[set-max] default = 500
MIN=$7      #[set-min] default = 5
SEED=$8     #[rnd_seed] number >= 149 digits

GSEA4="/home/samkyy/bin/GSEA_4.1.0"
cwd=$(pwd)

cd ${cwd}
echo ${cwd}

${GSEA4}/gsea-cli.sh GSEA -res ${RES} -cls ${PHENO} -gmx ${GENESET} -collapse No_Collapse -permute ${PERM} -nperm 1000 -scoring_scheme weighted -metric Signal2Noise -sort real -order descending -set_max ${MAX} -set_min ${MIN} -mode Max_probe -norm meandiv -rnd_type no_balance -include_only_symbols true -make_sets true -median false -num 200 -plot_top_x 20 -rnd_seed ${SEED} -save_rnd_lists false -zip_report false -rpt_label ${LABEL}

echo "GSEA run complete for ${LABEL}"

# Example: -res /home/samkyy/data/2021-05-06_gsea/GE_gsea_table.txt -cls /home/samkyy/data/2021-05-06_gsea/phenotypes.cls#cluster0_vs_REST -gmx /home/samkyy/data/2021-05-06_gsea/genesets.gmx -collapse false -permute gene_set -nperm 1000 -scoring_scheme weighted -metric Signal2Noise -sort real -order descending -set_max 500 -set-min 5 -mode Max_probe -norm meandiv -rnd_type no_balance -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -zip_report false -rpt_label cluster0_VarheekWang_GE -gui false
