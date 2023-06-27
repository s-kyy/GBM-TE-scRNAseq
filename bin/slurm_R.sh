#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

module load r/4.0.2
module load gcc/9.3.0
module load gdal/3.2.3

## test number of arguements passed

if [ $# -eq 0 ]
then
    echo "No arguments supplied"
elif [ $# -eq 1 ]
    DATAPATH=${1}
    Rscript ~/scratch/gete-gbm/bin/r_monocle3.R ${DATAPATH} &> ~/scratch/temp/outs
else
    DATAPATH=${1}
    OUTPUTPATH=${2}
    Rscript ~/scratch/gete-gbm/bin/r_monocle3.R ${DATAPATH} ${OUTPUTPATH}  &> ~/scratch/temp/outs
fi