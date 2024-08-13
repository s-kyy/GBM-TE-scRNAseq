#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --array=0-10
#SBATCH --mail-user=samantha.y.twentyfourteen@gmail.com
#SBATCH --mail-type=ALL

module load java/11.0.22

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster0_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster0_c7-immune_GTE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster1_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster1_c7-immune_GTE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster2_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster2_c7-immune_GTE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster3_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster3_c7-immune_GTE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 4 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster4_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster4_c7-immune_GTE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 5 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster5_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster5_c7-immune_GTE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 6 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster6_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster6_c7-immune_GTE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 7 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster7_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster7_c7-immune_GTE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 8 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster8_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster8_c7-immune_GTE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 9 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster9_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster9_c7-immune_GTE_table 500 5 160 
fi

if [[ $SLURM_ARRAY_TASK_ID == 10 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster10_versus_REST /home/samkyy/scratch/gsea/genesets/c7.immunesigdb.v2023.2.Hs.symbols.gmt gene_set cluster10_c7-immune_GTE_table 500 5 160 
fi

