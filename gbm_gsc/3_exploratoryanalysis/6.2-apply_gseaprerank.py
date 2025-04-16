#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-

# Python v3.11.5 

# Author: Samantha Yuen
# Adapted from: apply_cellranger6.py

# Project: RetroTransposon scRNA-seq analysis on adult human glioblastoma 

# Example
# python3 ~/scripts/6.2-apply_gseaprerank.py [ranked_genes_dir] [geneset_dir] [rpt_label] [set-max] [set-min] [rnd_seed] [out_dir]

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
tmpdir = os.path.join(cwdir,"tmp_scripts")
Path(tmpdir).mkdir(parents=False, exist_ok=True)
outdir = os.path.join(cwdir,"gsea_results")

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
    description='Produce script for GSEAPrerank on scRNA-seq clusters',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument('-i', nargs=1, help='[ranked_genes_dir]', type=parse_path, dest="rnk", required = True)
parser.add_argument('-g', nargs=1, help='[geneset_dir]', type=parse_path, dest="sets", required = True)
parser.add_argument('-l',  nargs=1, help='[rpt_label]', type=str, dest="label", required = True)
parser.add_argument('-mx', nargs=1, help='[set_max]', default=500, type=int, dest="max", required = False)
parser.add_argument('-mn', nargs=1, help='[set_min]', default=5, type=int, dest="min", required = False)
parser.add_argument('-r', nargs=1, help='[rnd_seed]', default=149, type=int, dest="seed", required = False)
parser.add_argument('-o', nargs=1, help='[out_dir]', default=outdir, type=parse_parent, dest="out", required = False)
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

rnk_dir = returnValue(args.rnk)
geneset_dir = returnValue(args.sets) 
outdir = returnValue(args.out) 
Path(outdir).mkdir(parents=False, exist_ok=True)

rpt_label = returnValue(args.label)
set_max = str(returnValue(args.max))
set_min = str(returnValue(args.min))
seed = str(returnValue(args.seed))

# ==============================================================================
# Extract list of filenames 
# ==============================================================================

def getFilePaths(dir, file_ext):
    """
    return a list of file paths with desired extension. 
    file paths are renamed if it contains spaces
    """
    file_list = []
    for root, dirs, filenames in os.walk(dir):
        for filename in filenames:
            if filename.endswith(file_ext):
                
                filename_renamed = filename.replace(" ","")
                filename_renamed = filename_renamed.replace("(","")
                filename_renamed = filename_renamed.replace(")","")
                if (filename_renamed != filename): 
                    os.rename(
                        os.path.join(root, filename),
                        os.path.join(root, filename_renamed)
                    )
                filepath = os.path.join(root, filename_renamed)

                file_list.append(filepath)
    return file_list

rnk_files = getFilePaths(rnk_dir, (".rnk"))
print(rnk_files)

geneset_files = getFilePaths(geneset_dir, (".gmx",".gmt"))
print(geneset_files)

# # ==============================================================================
# # Parse rnk file names to generate report labels
# # ==============================================================================

def split_str(seperator_list, str_to_split):
    # create regular expression dynamically
    regex = '|'.join(map(re.escape, seperator_list))
    return re.split(regex, str_to_split)

separators = "_"
sample_names = []
for rnk_file in rnk_files:
    seperated = split_str( separators, os.path.splitext(os.path.basename(rnk_file))[0] ) # split rnk_file name excluding file extension 
    sample_name = seperated[-2] + '_' + seperated[-1] # use last two elements as label name
    sample_names.append( sample_name.replace(" ", "") ) # strip whitespaces
print(sample_names)

# ==============================================================================
# Write GSEA script as slurm job
# ==============================================================================

# Example filename format: <cwd>/tmp_scripts/gbm_ge_gsea_YYYY-mm-dd_HhMm.sh
out_script = os.path.join(tmpdir, rpt_label + '_gsea_' + datetime.now().strftime("%Y-%m-%d_%Hh%Mm") + ".sh")
print(out_script)

num_jobs = (len(geneset_files)*len(rnk_files)) 

try:

    with open(out_script, 'w') as script:
        script.write("#!/bin/bash\n")
        print("Wrote shebang")
        script.write("#SBATCH --time=00:5:00\n")
        script.write("#SBATCH --account=xxx\n")
        script.write("#SBATCH --ntasks=1\n")
        script.write("#SBATCH --cpus-per-task=1\n")
        script.write("#SBATCH --mem-per-cpu=1G\n")
        script.write("#SBATCH --array=0-")
        script.write(str( num_jobs -1 ) + "\n")  
        print("Wrote script paramters")

except IOError:

    sys.exit('Cannot append in tmp SLURM file', 0)
    
# Loop through all combinations of expression datasets and genesets to create GSEA scripts. 

i = 0 

for geneset_file in geneset_files:

    sample_name_id = 0

    for rnk_file in rnk_files:

        try:
            
            with open(out_script, 'a') as script:

                ## Example Output ##
                # if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then ... fi

                print("writing job number ", str(i))

                script.write("if [[ $SLURM_ARRAY_TASK_ID == " + str(i) + " ]] ; then\n")
                script.write("~/scratch/runs/3_exploratoryanalysis/6.3-run_gseaprerank.sh \"" +
                                rnk_file + "\" \"" +    # RES=$1      #[expression_dataset] filepath
                                # "#" + sample_names[sample_name_id] + "_versus_REST " +
                                geneset_file + "\" " +     # GENESET=$2  #[geneset] filepath
                                rpt_label + "_" + sample_names[sample_name_id] + "_" + os.path.splitext(os.path.basename(geneset_file))[0] + " " +
                                                    # LABEL=$3    # Report Label
                                set_max + " " +     # MAX=$4      #[set-max] default = 500
                                set_min + " " +     # MIN=$5      #[set-min] default = 5
                                seed + " " +        # SEED=$6     #[rnd_seed] number >= 149 digits
                                outdir + " " +      # OUT=$7      #[output_dir] 
                                "\n")
                script.write("fi\n\n")

            i+=1 
            sample_name_id+=1
        
        except IOError:
        
            sys.exit('Cannot append in tmp SLURM file', 0)

print("Writing to script complete")