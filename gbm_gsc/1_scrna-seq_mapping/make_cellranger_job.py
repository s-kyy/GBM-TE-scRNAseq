#!/usr/bin/env python3.9 
# -*- coding: utf-8 -*- 


# Python v3.9.2 

# Author: Samantha Yuen
# Adapted from: apply_cellranger3.pl by Dr. Yoshiaki Tanaka
# Created: Thursday Mar 04, 2021
# Modified: Friday Nov 15, 2024

# Purpose: scRNA-seq read mapping and counting with cellranger from fastq files. 

# Example
# python3 ~/script/make_cellranger_job.py [fastqdir] [output_dir] [ref_genome] [localcores] [mempercore] [v3-v6] [output_name]

# Imports 
import sys, os, io, time, glob
from os import path
from datetime import datetime

# Global variables 
cwdir = os.getcwd()  # current working directory
tmpdir = "." # slurm script output directory
fastqdir = ""

### function: Verify paramater values
def verifyArgs():

    # Check for number of arguements
    if (not (len(sys.argv[1:]) == 7)):
        sys.exit("Usage: " + os.path.basename(sys.argv[0]) + " [fastqdir] [output_dir] [ref_genome] [localcores] [mempercore] [v3-v6] [output_name]\n")

    #Verify directories and values. 
    if (not os.path.exists(sys.argv[1])):
        sys.exit('Path to fastq files do not exist, terminating program...')
    if (not os.path.exists(sys.argv[2])):
        sys.exit('Path to output directory does not exist, terminating program...')
    if (not os.path.exists(sys.argv[3])):
        sys.exit('Path to reference genome does not exist, terminating program...')
    try:
        int(sys.argv[4])
    except ValueError:
        sys.exit('Invalid localcores size, must be an int, terminating program...')
    try:
        int(sys.argv[5])
    except ValueError:
        sys.exit('Invalid mempercore size, must be an int, terminating program...')
    try: 
        int(sys.argv[6])

        if (int(sys.argv[6]) != 3 and int(sys.argv[6]) != 6):
            sys.exit('Invalid v3-v6 option, value must be 3 or 6, terminating program...')
    except ValueError:
        sys.exit('Invalid v3-v6 option, must be and int, terminating program...')

def findFolders(fastqdir):

    folderPaths = [] # list for fastq directories

    # print("in findFolders def for:", fastqdir)

    # for filePath in glob.glob(fastqdir + "/data_*/*/*_I1_*.fastq.gz"):
        #glob a list of all paths that contain *_I1_*.fastq.gz files
        #ex: fastqdir = RetroTransposonAnalysis/data00/bamtofastqfolder/*_I1_*.fastq.gz 

    for filePath in glob.glob(fastqdir + "/SRR*/*_R1_001.fastq.gz"):
        # ex: fastqdir = gete-gbm/data/SRR10353960/GBM27_*_R1_001.fastq.gz

        print ("a filePath was found:", os.path.dirname(filePath), "\n")
    
        folderPaths.append(os.path.dirname(filePath)) #strip the name of the files
        
        # print("Length of folderPaths list:", len(folderPaths))

        if not os.path.isdir(os.path.dirname(filePath)):
            print("directory does not exist\n")
            sys.exit()
    
    result = [*{*folderPaths}]
    return  result #return unique list of folders

### Launcher 
def main(): 
    # Arguements passed from command line
    fastqdir = sys.argv[1]  # directory of fastq files
    outdir = sys.argv[2]    # output of map cell/genes count matrix (named map/)
    ref_gen = sys.argv[3]   # reference genome
    localcores = sys.argv[4]  # local cores value used in SBATCH (e.g. 8, 12, 16)
                            # larger the value the higher greater the faster the process 
    mempercore = sys.argv[5] # mempercore value used in SBATCH (e.g. 10, 12)
    version = int(sys.argv[6])   # version of cellranger to run (3 for retrotransposons referencen genes, else 6)
    outName = sys.argv[7]
    
    print(fastqdir) # folder containing fastq folders we are processing. 


    print("finding folders:\n")

    os.chdir(fastqdir) #change directories containing the fastq files

    uniqueFolderPaths = findFolders(fastqdir)
    print("uniqueFolderPaths size", len(uniqueFolderPaths))

    # move to workdirectory for slurm output script
    if (not os.path.exists(tmpdir)):
        try: 
            os.mkdir(tmpdir)
        except OSError:
            print("Error making tmpdir for slurm output/")
            sys.exit()
    os.chdir(tmpdir)

    # make tmp file for slurm bash script
    tmp = tmpdir + "/" + os.path.basename(outName) +  datetime.now().strftime("%Y-%m-%d_%Hh%Mm") + ".sh"
    print("Output slurm file:", os.path.basename(tmp))

    try: 
        with open(tmp, 'w') as slurm:
            slurm.write("#!/bin/bash\n")
            slurm.write("#SBATCH --time=36:00:00\n")
            slurm.write("#SBATCH --account=xxx\n")
            slurm.write("#SBATCH --ntasks=1\n")
            slurm.write("#SBATCH --cpus-per-task=" + localcores + "\n")
            slurm.write("#SBATCH --mem-per-cpu=" + mempercore + "G\n")
            # slurm.write("#SBATCH --mail-user=xxx\n")
            # slurm.write("#SBATCH --mail-type=ALL\n")
            slurm.write("#SBATCH --array=0-")
            slurm.write(str(len(uniqueFolderPaths)-1) + "\n")

            slurm.write("cd " + outdir + "\n") #should be cd-ing into ./script, but since I put absolute path to run_cellranger6.sh, it's fine
            
            # array output: if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then ... fi
            # ${CELLRANGER}/cellranger count ${OPTION} --localcores=${CPU} --id=${OUTPUT} --transcriptome=${REF} --fastq=${INPUT} 
            for i in range(len(uniqueFolderPaths)): # range: 0-n
                slurm.write("if [[ $SLURM_ARRAY_TASK_ID == " + str(i) + " ]] ; then\n")
                if (version == 3):
                    slurm.write("../run_cellranger3.sh " + 
                                uniqueFolderPaths[i] + " " +    # --fastq=${INPUT}
                                "map_" + os.path.basename(uniqueFolderPaths[i]) + " " +     
                                                                # --id=${OUTPUT} output directory
                                outdir + " " +                  # cd into containing all output directories in $DIR
                                ref_gen + " " +                 # --transcriptome=${REF}
                                localcores + " " +              # --localcores=${CPU}
                                mempercore + " " +              # --mempercore=${OPTION}
                                "\n")
                else:
                    slurm.write("../run_cellranger6.sh " + 
                                    uniqueFolderPaths[i] + " " +    # --fastq=${INPUT}
                                    "map_" + os.path.basename(uniqueFolderPaths[i]) + " " +     # --id=${OUTPUT} output directory
                                    outdir + " " +                # cd into containing all output directories in $DIR
                                    ref_gen + " " +                 # --transcriptome=${REF}
                                    localcores + " " +              # --localcores=${CPU}
                                    mempercore + " " +              # --mempercore=${OPTION}
                                    "\n")

                slurm.write("fi\n\n")
    except IOError:
        print('Cannot append in tmp SLURM file')    
        sys.exit()

### main
if __name__ == "__main__": 

    verifyArgs()
    main() 
