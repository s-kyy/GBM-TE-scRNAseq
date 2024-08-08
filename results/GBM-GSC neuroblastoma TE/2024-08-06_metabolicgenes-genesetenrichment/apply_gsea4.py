#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-

# Python v3.11.5 

# Author: Samantha Yuen
# Adapted from: apply_cellranger6.py
# Created: Tuesday May 18, 2021
# Modified: Tuesday August 6, 2024

# Project: RetroTransposon scRNA-seq analysis on adult human glioblastoma 

# Example
# python3 ~/scripts/apply_gsea4.py [expression_dataset] [phenotype_label] [geneset] [rpt_label] [permute] [set-max] [set-min] [rnd_seed]
# python3 ~/scripts/apply_gsea4.py /home/samkyy/data/2021-05-06_gsea/GE_gsea_table.txt /home/samkyy/data/2021-05-06_gsea/phenotypes.cls /home/samkyy/data/2021-05-06_gsea/genesets.gmx VarheekWang gene_set 500 5 160

# Imports 
import sys, os, io, time, re
from os import path
from datetime import datetime

# Global variables
cwdir = os.getcwd()  # current working directory
tmpdir = "/home/samkyy/scratch/slurm_tmp" # slurm script output directory
separators = "# ", "\n", " "

### function: Verify paramater values
def verifyArgs():

    # Check for number of arguements
    if (not (len(sys.argv[1:]) >= 4)):
        sys.exit("Usage: " + os.path.basename(sys.argv[0]) + " [expression_dataset] [phenotype_label] [geneset] [rpt_label] [permute] [set-max] [set-min] [rnd_seed]\n")

    #Verify directories and values. 
    if (not os.path.exists(sys.argv[1])):
        sys.exit('Path to EXPRESSION DATASET does not exist, terminating program...')
    if (not os.path.exists(sys.argv[2])):
        sys.exit('Path to PHENOTYPE LABELS does not exist, terminating program...')
    if (not os.path.exists(sys.argv[3])):
        sys.exit('Path to GENESET does not exist, terminating program...')
    try:
        str(sys.argv[4])
    except ValueError:
        sys.exit('Invalid rpt_label value, must be a string, terminating program...')
    try:
        str(sys.argv[5])
    except ValueError:
        sys.exit('Invalid permute value, must be "gene_set" or "phenotype", terminating program...')
    try:
        int(sys.argv[6])
    except ValueError:
        sys.exit('Invalid set-max value, must be an int, terminating program...')
    try:
        int(sys.argv[6])
    except ValueError:
        sys.exit('Invalid set-min value, must be an int, terminating program...')
    try:
        int(sys.argv[7])
    except ValueError:
        sys.exit('Invalid rnd_seed value, must be an int, terminating program...')
    
def read_file(fname):
    with open(fname, "r") as file:
        data = file.readlines()
    return data

def split_str(sepr_list, str_to_split):
    # create regular expression dynamically
    regular_exp = '|'.join(map(re.escape, sepr_list))
    return re.split(regular_exp, str_to_split)

### Launcher
def main():
    # Arguements passed from command line
    res = sys.argv[1]
    pheno = sys.argv[2]
    geneset = sys.argv[3]
    rpt_label = sys.argv[4]
    perm = sys.argv[5]
    set_max = sys.argv[6]
    set_min = sys.argv[7]
    seed = sys.argv[8] if int(sys.argv[8]) >= 149 else str(149)

    # Extract list of sample names from .cls phenotype lables
    pheno_raw = read_file(pheno)
    snames = [x  for x in split_str(separators,pheno_raw[1]) if x] # generate a list that consists of x for every element in strings if x actually contains something
    print(snames)

    # make tmp file for slurm bash script
    tmp = tmpdir + "/" + os.path.splitext(os.path.basename(res))[0] + rpt_label + datetime.now().strftime("%Y-%m-%d_%Hh%Mm") + ".sh"
        # /home/samkyy/scratch/slurm_tmp/GE_gsea_tableYYYY-mm-dd_HhMm
    print("Output slurm file:", os.path.basename(tmp))

    try:
        with open(tmp, 'w') as slurm:
            slurm.write("#!/bin/bash\n")
            print("Wrote shebang")
            slurm.write("#SBATCH --time=00:05:00\n")
            slurm.write("#SBATCH --account=def-ytanaka\n")
            slurm.write("#SBATCH --ntasks=1\n")
            slurm.write("#SBATCH --cpus-per-task=1\n")
            slurm.write("#SBATCH --mem-per-cpu=2G\n")
            slurm.write("#SBATCH --array=0-")
            slurm.write(str(len(snames)-1) + "\n")
            slurm.write("#SBATCH --mail-user=samantha.y.twentyfourteen@gmail.com\n")
            slurm.write("#SBATCH --mail-type=ALL\n")
            print("Wrote slurm paramters")

            ## Example Output ##
            # if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then ... fi

            for i in range(len(snames)): # range: 1-22
                print("writing array element number ", str(i))
                slurm.write("if [[ $SLURM_ARRAY_TASK_ID == " + str(i) + " ]] ; then\n")
                slurm.write("/home/samkyy/scratch/gete-gbm/bin/run_gsea4.sh " +
                                res + " " +         # RES=$1      #[expression_dataset] filepath
                                pheno +             # PHENO=$2    #[phenotype_label] filepath
                                "#" + snames[i] + "_versus_REST " +
                                                    # EXP=$5      # Samples to compare (Include the pound symbol) ex: sample_vs_REST
                                geneset + " " +     # GENESET=$3  #[geneset] filepath
                                perm + " " +        # PERM=$4     #[permute] gene-set or phenotype
                                snames[i] + "_" + rpt_label + "_" + os.path.splitext(os.path.basename(res))[0] + " " +
                                                    # LABEL=$6    # Report Label
                                set_max + " " +     # MAX=$7      #[set-max] default = 500
                                set_min + " " +     # MIN=$8      #[set-min] default = 5
                                seed + " " +        # SEED=$9     #[rnd_seed] number >= 149 digits
                                "\n")
                slurm.write("fi\n\n")

            print("Writing to slurm file complete")
    except IOError:
        sys.exit('Cannot append in tmp SLURM file', 0)

### main
if __name__== "__main__":

    startTime = time.time() #float

    verifyArgs()
    main()

    endTime = time.time() #float

    print("start:", startTime, "\n")
    print("end:", endTime, "\n")
    print("taken:", endTime-startTime, "\n")
