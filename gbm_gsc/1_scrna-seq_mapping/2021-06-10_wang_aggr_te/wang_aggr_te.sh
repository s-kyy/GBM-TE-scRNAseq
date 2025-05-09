#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=15G

scriptdir="$(dirname "$0")"
cd "${scriptdir}/2021-06-10_wang_aggr_te"
echo "$scriptdir" 

~/bin/cellranger-3.0.2/cellranger aggr --id=aggr --csv=./aggr_te.csv
_ge.csv --normalize=none --nosecondary