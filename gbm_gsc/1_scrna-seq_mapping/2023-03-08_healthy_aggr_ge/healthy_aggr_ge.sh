#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G

scriptdir="$(dirname "$0")"
cd "${scriptdir}/2023-03-08_healthy_aggr_ge"
echo "$scriptdir" 

~/bin/cellranger-6.0.0/cellranger aggr --id=aggr --csv=./aggr_ge.csv
_ge.csv --normalize=none --nosecondary