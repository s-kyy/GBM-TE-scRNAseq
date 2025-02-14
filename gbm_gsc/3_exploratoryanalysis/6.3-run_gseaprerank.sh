#!/bin/sh

RNK=$1      #[expression_dataset] filepath
GENESET=$2  #[geneset] filepath

LABEL=$3    # Report Label

MAX=$4      #[set-max] default = 500
MIN=$5      #[set-min] default = 5
SEED=$6     #[rnd_seed] number >= 149 digits
OUT=$7      #[out] output directory

module load java/17.0.2

GSEA4="/home/samkyy/bin/GSEA_4.1.0"
cwd=$(pwd)

cd ${cwd}
echo ${cwd}

${GSEA4}/gsea-cli.sh GSEAPreranked -rnk ${RNK} -gmx ${GENESET} -collapse No_Collapse -nperm 1000 -scoring_scheme weighted -set_max ${MAX} -set_min ${MIN} -mode Max_probe -norm meandiv -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed ${SEED} -zip_report false -rpt_label ${LABEL} -create_svgs false -out ${OUT}

echo "GSEAPreranked run complete for ${RNK} as ${LABEL} in ${OUT}"
