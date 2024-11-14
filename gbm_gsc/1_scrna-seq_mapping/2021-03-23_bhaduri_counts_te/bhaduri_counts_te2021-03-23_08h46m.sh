#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G
#SBATCH --array=0-11
cd /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs
if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_06/SF11232_1_MissingLibrary_1_HKMVGBCXY mapRT_0 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_05/SF22215_GBM2_MissingLibrary_1_CAUFRANXX mapRT_1 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_01/SF11159_2_MissingLibrary_1_HFMLTBCXY mapRT_2 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_8/Tumor2_1_MissingLibrary_1_CB64GANXX mapRT_3 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 4 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_07/SF11232_2_MissingLibrary_1_HKMVGBCXY mapRT_4 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 5 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_02/SF11209_GBM1_MissingLibrary_1_CAUFRANXX mapRT_5 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 6 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_03/SF11209_GBM2_MissingLibrary_1_CAUFRANXX mapRT_6 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 7 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_00/SF11159_1_MissingLibrary_1_HFMLTBCXY mapRT_7 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 8 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_04/SF22215_GBM1_MissingLibrary_1_CAUFRANXX mapRT_8 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 9 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_10/TQ1_MissingLibrary_1_CAT6PANXX mapRT_9 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 10 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_11/TQ2_MissingLibrary_1_CAT6PANXX mapRT_10 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

if [[ $SLURM_ARRAY_TASK_ID == 11 ]] ; then
/home/samkyy/bin/run_cellranger3.sh /home/samkyy/scratch/RetroTransposonAnalysis/data_09/Tumor2_2_MissingLibrary_1_CB64GANXX mapRT_11 /home/samkyy/scratch/RetroTransposonAnalysis/cellranger_outputs /home/samkyy/projects/def-ytanaka/common/refdata_GRCh38-TE_v2 8 15 
fi

