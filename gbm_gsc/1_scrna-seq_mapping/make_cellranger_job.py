#!/usr/bin/env python3.9 
# -*- coding: utf-8 -*- 


# Python v3.9.2 

# Author: Samantha Yuen
# Adapted from: apply_cellranger3.pl by Dr. Yoshiaki Tanaka
# Created: Thursday Mar 04, 2021
# Modified: Friday Nov 15, 2024

# Purpose: scRNA-seq read mapping and counting with cellranger from fastq files. 

# Example
# python3 ~/script/make_cellranger_job.py [fastqdir] [output_dir] [ref_genome] [localcores] [mempercore] [v3-v6]

# ==============================================================================
# Import Librairies
# ==============================================================================
import sys, os, glob, argparse
from datetime import datetime
# ==============================================================================
# Global variables 
# ==============================================================================
cwdir = os.getcwd()  # current working directory

# ==============================================================================
# Verify arguments
# ==============================================================================
def parse_path(value):
    if os.path.exists(value):
        return value
    else:
        sys.exit('Paths do not exist. Terminating program...')

def returnValue(obj):
    """
    return the first item in an argparse variable for optional arguments 
    """
    if isinstance(obj, list):
        return obj[0]
    else:
        return obj

parser = argparse.ArgumentParser(description='Produce script for cellranger count function',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i',nargs=1,help='fastq filepath',type=parse_path,dest="input_folder",required=True)
parser.add_argument('-o',nargs=1,help='output filepath for cellranger count (ie. "./")',type=parse_path,dest="out_folder",required=True) 
parser.add_argument('-r',nargs=1,help='reference annotations filepath', type=parse_path, dest="ref", required=True)
parser.add_argument('-c', nargs=1, default=8, help='local cores value used in SBATCH (e.g. 8, 12, 16). larger the value the higher greater the faster the process ', type=int, dest="cores")
parser.add_argument('-m', nargs=1, default=10, help="mempercore value used in SBATCH (e.g. 10, 12, 15)", type=int, dest="mem")
parser.add_argument('-t', nargs=1, default=6, help="cellranger version to use. (6 or 3)", type=int, dest="type")
parser.add_argument('-p', nargs=1, default=0, help="SF... = 0 (default), SRR... = 1", type=int, dest="pattern")

args = parser.parse_args()

# ==============================================================================
# Cleanup args values
# ==============================================================================
fastqdir = returnValue(args.input_folder)
outdir = returnValue(args.out_folder) 
ref_gen = returnValue(args.ref) 
localcores = returnValue(args.cores)
mempercore = returnValue(args.mem)
cellranger_type = returnValue(args.type)
pattern = returnValue(args.pattern)
print(fastqdir) 

# ==============================================================================
# Extract list of unique folders
# ==============================================================================
folderPaths = [] # list for fastq directories paths
foldernames = [] # list for fastq directory names
sample_name_pattern = ""

#change to directory containing the fastq files
os.chdir(fastqdir) 

print('Extracting file paths by sample')

# adjust globbed pattern based on the the fastq file pattern to the sample prefix
if pattern == 0 :
    sample_name_pattern = "/SF*_*/*/*_I1_*.fastq.gz"
    #ex: fastqdir = fastqdir/samplename/bamtofastqfolder/*_I1_*.fastq.gz 
else:
    sample_name_pattern = "/SRR*/*_R1_001_.fastq.gz"

for filePath in glob.glob(fastqdir + sample_name_pattern): 
    #ex: fastqdir = fastqdir/samplename/bamtofastqfolder/*_I1_*.fastq.gz 
    print (os.path.dirname(filePath), "\n")
    folderPaths.append(os.path.dirname(filePath)) #strip the name of the files
    foldernames.append(os.path.basename(os.path.dirname(os.path.dirname(filePath))))

uniqueFolderPaths = sorted([*{*folderPaths}])
foldernames = sorted([*{*foldernames}])
print(foldernames)
print("uniqueFolderPaths size", len(uniqueFolderPaths))

# ==============================================================================
# Write script
# ==============================================================================

# move to workdirectory for slurm output script
os.chdir(cwdir) 

# make tmp file for slurm bash script
tmp = os.path.join(cwdir, os.path.basename(fastqdir), datetime.now().strftime("%Y-%m-%d_%Hh%Mm"), ".sh")
    # ./<fastqdir>/YYYY-mm-dd_HhMm.sh
print("Output slurm file:", tmp)

try: 
    with open(tmp, 'w') as slurm:
        slurm.write("#!/bin/bash\n")
        slurm.write("#SBATCH --time=36:00:00\n")
        slurm.write("#SBATCH --account=xxx\n")
        slurm.write("#SBATCH --ntasks=1\n")
        slurm.write("#SBATCH --cpus-per-task=" + str(localcores) + "\n")
        slurm.write("#SBATCH --mem-per-cpu=" + str(mempercore) + "G\n")
        slurm.write("#SBATCH --array=0-")
        slurm.write(str(len(uniqueFolderPaths)-1) + "\n")
        # slurm.write("#SBATCH --mail-user=\n")
        # slurm.write("#SBATCH --mail-type=ALL\n")

        slurm.write("cd " + outdir + "\n") #should be cd-ing into ./script, but since I put absolute path to run_cellranger6.sh, it's fine
        
        # Output Example:
        # array output: if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then 
        # ${CELLRANGER}/cellranger count ${OPTION} --localcores=${CPU} --id=${OUTPUT} --transcriptome=${REF} --fastq=${INPUT} 
        # fi
        
        for i in range(len(uniqueFolderPaths)): # range: 0-5
            slurm.write("if [[ $SLURM_ARRAY_TASK_ID == " + str(i) + " ]] ; then\n")
            slurm.write("../run_cellranger" + str(cellranger_type) + ".sh " + 
                            uniqueFolderPaths[i] + " " +    # --fastq=${INPUT}
                            "map_" + foldernames[i] + " " +         # --id=${OUTPUT} output directory
                            outdir + " " +                  # cd into containing all output directories in $DIR
                            ref_gen + " " +                 # --transcriptome=${REF}
                            str(localcores) + " " +              # --localcores=${CPU}
                            str(mempercore) + " " +              # --mempercore=${OPTION}
                            "\n")

            slurm.write("fi\n\n")
except IOError:
    print('Cannot append in tmp SLURM file')    
    sys.exit()