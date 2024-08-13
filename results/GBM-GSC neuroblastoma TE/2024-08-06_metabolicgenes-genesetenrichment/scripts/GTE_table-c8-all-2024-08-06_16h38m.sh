#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --array=0
#SBATCH --mail-user=samantha.y.twentyfourteen@gmail.com
#SBATCH --mail-type=ALL

module load java/11.0.22

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
/home/samkyy/scratch/gsea/run_gsea4.sh /home/samkyy/scratch/gsea/GTE_table.gct /home/samkyy/scratch/gsea/phenotypes_gte.cls#cluster0_versus_REST /home/samkyy/scratch/gsea/genesets/c8.all.v2023.2.Hs.symbols.gmt gene_set cluster0_c8-all_GTE_table 500 5 160 
fi

