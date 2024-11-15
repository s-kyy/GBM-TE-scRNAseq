#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

# Python v3.9.2 

# Author: Samantha Yuen
# Modified: Tuesday Mar 07, 2023     

# ==============================================================================
# Import Librairies
# ==============================================================================
import sys, os, glob, argparse
from os import path
from datetime import datetime

# ==============================================================================
# Global variables 
# ==============================================================================
cwdir = os.getcwd()  # current working directory

# ==============================================================================
# Verify argument
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
parser.add_argument('-c', nargs=1, default=8, help='local cores value used in SBATCH (e.g. 8, 12, 16). larger the value the higher greater the faster the process ', type=int, dest="cores")
parser.add_argument('-m', nargs=1, default=10, help="mempercore value used in SBATCH (e.g. 10, 12, 15)", type=int, dest="mem")

args = parser.parse_args()

# ==============================================================================
# Cleanup args values
# ==============================================================================
fastqdir = returnValue(args.input_folder)
localcores = returnValue(args.cores)
mempercore = returnValue(args.mem)
print(fastqdir) 

# ==============================================================================
# Extract list of unique folders
# ==============================================================================
folderPaths = [] # list for fastq directories paths

os.chdir(fastqdir) #change to directory containing the fastq files

print('Extracting file paths by sample')
for filePath in glob.glob(fastqdir + "/SRR*/*_R*.fastq"):
    print (filePath, "\n")
    folderPaths.append(filePath) #strip the name of the files

uniqueFolderPaths = sorted([*{*folderPaths}])
print("uniqueFolderPaths size", len(uniqueFolderPaths))
print("uniqueFolderPaths", uniqueFolderPaths)

# move to workdirectory for slurm output script
os.chdir(cwdir) #change to directory containing the fastq files

# make tmp file for slurm bash script
tmp = cwdir + "/pigz_" + datetime.now().strftime("%Y-%m-%d_%Hh%Mm") + ".sh"
    # /home/samkyy/scratch/slurm_tmp/fastqdir/YYYY-mm-dd_HhMm
print("Output slurm file:", os.path.basename(tmp))

try: 
    with open(tmp, 'w') as slurm:
        slurm.write("#!/bin/bash\n")
        slurm.write("#SBATCH --time=2:00:00\n")
        slurm.write("#SBATCH --account=xxx\n")
        slurm.write("#SBATCH --ntasks=1\n")
        slurm.write("#SBATCH --cpus-per-task=" + str(localcores) + "\n")
        slurm.write("#SBATCH --mem-per-cpu=" + str(mempercore) + "G\n")
        slurm.write("#SBATCH --array=0-")
        slurm.write(str(len(uniqueFolderPaths)-1) + "\n")
        slurm.write("#SBATCH --mail-user=\n")
        slurm.write("#SBATCH --mail-type=ALL\n")

        slurm.write("cd " + fastqdir + "\n")
        
        for i in range(len(uniqueFolderPaths)): # range: 0-5
            slurm.write("if [[ $SLURM_ARRAY_TASK_ID == " + str(i) + " ]] ; then\n")
            slurm.write("pigz -8 -b 1024 -vk -p " + 
                            str(localcores) + " " +
                            uniqueFolderPaths[i] + " " +
                            "\n")
            slurm.write("pigz -vt " +
                            uniqueFolderPaths[i] + ".gz" + 
                            "\n")
            slurm.write("fi\n\n")
except IOError:
    print('Cannot append in tmp SLURM file')    
    sys.exit()