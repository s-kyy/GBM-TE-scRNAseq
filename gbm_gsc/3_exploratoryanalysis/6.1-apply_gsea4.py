#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-

# Python v3.11.5 

# Author: Samantha Yuen
# Adapted from: apply_cellranger6.py
# Created: Tuesday May 18, 2021
# Modified: Tuesday Feb 3, 2024

# Project: RetroTransposon scRNA-seq analysis on adult human glioblastoma 

# Example
# python3 ~/scripts/apply_gsea4.py [expression_dataset] [phenotype_label] [geneset] [rpt_label] [set-max] [set-min] [rnd_seed]
# python3 ~/scripts/apply_gsea4.py ~/data/2021-05-06_gsea/GE_gsea_table.txt ~/data/2021-05-06_gsea/phenotypes.cls ~/data/2021-05-06_gsea/genesets.gmx VarheekWang 500 5 160

# ==============================================================================
# Imports 
# ==============================================================================
import sys, os, re, argparse
from os import path
from datetime import datetime
from pathlib import Path

# ==============================================================================
# Global variables 
# ==============================================================================
cwdir = os.getcwd()  # current working directory
tmpdir = Path(os.path.join(cwdir,"tmp_scripts")).mkdir(parents=False, exist_ok=True)
outdir = os.path.join(cwdir,"gsea_results")
separators = "# ", "\n", " "

# ==============================================================================
# Verify arguments
# ==============================================================================
def parse_path(value):
    """
    check if path exists, if not exit script 
    """
    if os.path.exists(value):
        return value
    else:
        sys.exit('Path ',value,' does not exist. Terminating program...')

def parse_parent(value):
    """
    check if path to parent directory exists, if not exit script 
    """
    if os.path.exists(path.dirname(value)):
        return value
    else:
        sys.exit('Path ',value,' does not exist. Terminating program...')

parser = argparse.ArgumentParser(
    description='Produce script for GSEA on scRNA-seq clusters',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument('-i', nargs=1, help='[expression_dataset_dir]', type=parse_path, dest="res", required = True)
parser.add_argument('-ph', nargs=1, help='[phenotype_label]', type=parse_path, dest="pheno", required = True)
parser.add_argument('-g', nargs=1, help='[geneset_dir]', type=parse_path, dest="sets", required = True)
parser.add_argument('-l',  nargs=1, help='[rpt_label]', type=str, dest="label", required = True)
parser.add_argument('-mx', nargs=1, help='[set_max]', default=500, type=int, dest="max", required = True)
parser.add_argument('-mn', nargs=1, help='[set_min]', default=5, type=int, dest="min", required = True)
parser.add_argument('-r', nargs=1, help='[rnd_seed]', default=149, type=int, dest="seed", required = True)
parser.add_argument('-o', nargs=1, help='[out_dir]', default=outdir, type=parse_parent, dest="out", required = True)
args = parser.parse_args()

# ==============================================================================
# Cleanup args values
# ==============================================================================
def returnValue(obj):
    """
    return the first item in an argparse variable for optional arguments 
    """
    if isinstance(obj, list):
        return obj[0]
    else:
        return obj

res_dir = returnValue(args.res)
pheno = returnValue(args.pheno) 
geneset_dir = returnValue(args.sets) 
outdir = returnValue(args.out) 
Path(outdir).mkdir(parents=False, exist_ok=True)

rpt_label = returnValue(args.label)
set_max = returnValue(args.max)
set_min = returnValue(args.min)
seed = returnValue(args.seed)
perm = "gene_set"

# ==============================================================================
# Extract list of filenames from res_dir (.gct)
# ==============================================================================

res_files = []
for root, dirs, files in os.walk(res_dir):
    for file in files:
        if file.endswith(".gct"):
            res_files.append(file)
print(res_files)

# ==============================================================================
# Extract list of filenames from geneset_dir (.gmx, .gmt)
# ==============================================================================

geneset_files = []
for root, dirs, files in os.walk(geneset_dir):
    for file in files:
        if file.endswith(".gmx") or file.endswith(".gmt"):
            geneset_files.append(file)
print(geneset_files)

# ==============================================================================
# Extract list of sample names from .cls phenotype lables
# ==============================================================================

def read_file(file_name):
    with open(file_name, "r") as file:
        data = file.readlines()
    return data

def split_str(seperator_list, str_to_split):
    # create regular expression dynamically
    regex = '|'.join(map(re.escape, seperator_list))
    return re.split(regex, str_to_split)

pheno_raw = read_file(pheno)

# generate a list that consists of x for every element in strings if x actually contains something
sample_names = [x  for x in split_str(separators,pheno_raw[1]) if x] 
print(sample_names, " the first variable will be used as the analyzed sample")

# ==============================================================================
# Write GSEA script as slurm job
# ==============================================================================

# Example filename format: ~/scratch/tmp_scripts/gbm_ge_matrix_gsea_3ca_h_gbmsubtypes_YYYY-mm-dd_HhMm.sh
out_script = os.path.join(tmpdir, rpt_label + '_gsea_' + datetime.now().strftime("%Y-%m-%d_%Hh%Mm") + ".sh")
print(out_script)

sample_names = [value for value in sample_names if value !="rest"]

num_jobs = (len(geneset_files)*len(res_files)) # = number of loops and array size for slurm script

try:

    with open(out_script, 'w') as script:
        script.write("#!/bin/bash\n")
        print("Wrote shebang")
        script.write("#SBATCH --time=00:30:00\n")
        script.write("#SBATCH --account=xxx\n")
        script.write("#SBATCH --ntasks=1\n")
        script.write("#SBATCH --cpus-per-task=1\n")
        script.write("#SBATCH --mem-per-cpu=2G\n")
        script.write("#SBATCH --array=0-")
        script.write(str( num_jobs -1 ) + "\n")  
        print("Wrote script paramters")

except IOError:

    sys.exit('Cannot append in tmp SLURM file', 0)
    
# Loop through all combinations of expression datasets and genesets to create GSEA scripts. 

i = 0

for geneset_file in geneset_files:
    
    for res_file in res_files:

        try:
            
            with open(out_script, 'a') as script:

                ## Example Output ##
                # if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then ... fi

                print("writing array element number ", str(i))

                script.write("if [[ $SLURM_ARRAY_TASK_ID == " + str(i) + " ]] ; then\n")
                script.write("~/scratch/runs/3-exploratoryanaysis/6.2-run_gsea4.sh " +
                                res_file + " " +    # RES=$1      #[expression_dataset] filepath
                                pheno +             # PHENO=$2    #[phenotype_label] filepath
                                "#" + sample_names[0] + "_versus_REST " +
                                geneset_file + " " +     # GENESET=$3  #[geneset] filepath
                                perm + " " +        # PERM=$4     #[permute] gene-set or phenotype (GSEAPreranked is always by gene-set)
                                rpt_label + "_" + os.path.splitext(os.path.basename(res_file))[0] + " " +
                                                    # LABEL=$6    # Report Label
                                set_max + " " +     # MAX=$7      #[set-max] default = 500
                                set_min + " " +     # MIN=$8      #[set-min] default = 5
                                seed + " " +        # SEED=$9     #[rnd_seed] number >= 149 digits
                                outdir + " " +      # OUT=$9      #[output_dir] 
                                "\n")
                script.write("fi\n\n")

            i+=1 # increment job number
        
        except IOError:
        
            sys.exit('Cannot append in tmp SLURM file', 0)

print("Writing to script complete")