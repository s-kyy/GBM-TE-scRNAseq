#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --array=0-11
#SBATCH --mail-user=samantha.y.twentyfourteen@gmail.com
#SBATCH --mail-type=ALL
if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster0_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster0_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster1_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster1_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster2_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster2_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster3_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster3_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 4 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster4_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster4_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 5 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster5_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster5_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 6 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster6_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster6_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 7 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster7_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster7_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 8 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster8_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster8_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 9 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster9_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster9_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 10 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster10_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster10_c6-all_GE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 11 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GE_table.gct /home/samkyy/scratch/gsea/phenotypes_ge.cls#cluster11_versus_REST /home/samkyy/scratch/gsea/genesets/c6.all.v2023.2.Hs.symbols.gmt gene_set cluster11_c6-all_GE_table 500 5 160 
fi

