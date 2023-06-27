## Purpose

On 2023-02-28 Dr. Tanaka wanted to see the positively and negatively expressed genes in all clusters across both ge and ge+te annotated datasets. 

I looked at when the last time I ran this type of analysis, and realized I only ran it for the 2021-06-21 gbmgsc dataset at PC20r0.7 (`neuroblastoma.ipynb`, & `"D:\backup files\2021-11-11\Cluster\scratch\gete-gbm\results\GBM-GSC neuroblastoma\2021-06-21"`). Plus this only contained DEGs that were positively expressed. So, I ran this analysis again on 2023-02-28. 

## Methods + Notes

Ran sbatch with `R CMD BATCH` of two files to output csv.

- `./bin/r_gbmsc_findallmarkers.sh` --> `sbatch r_gbmsc_findallmarkers.sh`
- `./resultszxx/GBMGSCTE/r_findallmarkers_ge.sh`	
- `./results/GBMGSCTE/r_findallmarkers_gte.sh`
- Runs `FindAllMarkers()` function from [[Seurat]] and outputs the results as a csv file. 
- Job ID: `61400271_[0-1]` for `sbatch ./bin/r_gbmsc_findallmarkers.sh` -- Ran out of memory so I bumped it up to 16G in the next run
- Job ID: `61403098_[0-1]` for `sbatch ./bin/r_gbmsc_findallmarkers.sh` 
  - GE markers failed due to memory spike
  - GTE markers completed, ~~but returned an empty dataframe because none of the genes passed logfc.threshold > 0.25~~. The file wasn't empty, but there seems to be a small error that says makes it seem so based on the .Rout file.== If there are any issues from Dr. Tanaka I'll rerun it again with a lower fc threshold. 
- Job ID `61416756` ran `srun --time=4:00:0 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=32G --account=def-ytanaka R CMD BATCH ~/scratch/gete-gbm/results/GBMGSCTE/r_findallmarkers_ge.R` 
  - I only need less than an hour to run this command, but for good measure I'll put the full 2 hours in case it needs more time. 
  - ==GE markers script, Completed and was successful.==