#!usr/bin/env Rscript

# Saves output as a new directory inside the directory the script is run from.

.libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 2 filepath must be supplied: [xxx.rds] [cluster_col] [fig_dir_name]", call.=FALSE)
} else if (length(args)>=1) {
  
  # verify filepaths
  if (file.exists(args[1])) { 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    # parent_dir_name_obj <- basename(parent_dir_path_obj)
  } else {
    stop("Filepaths provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>1) {
    cluster_col <- args[2]
  } else {
    stop("Cluster Column not provided. Exiting...", call.=FALSE)
  }
  if (length(args)>2) {
    figs_dir_name <- args[3]
  } else {
    figs_dir_name <- "figs_markers"
  }

} else {
    stop("Error in calling R or filepaths. Exiting...", call.=FALSE)
}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
set.seed(108)

library(spatstat.core)
library(Seurat)
library(Matrix)
library(tidyverse) 
library(dplyr)
library(future)
library(future.apply)

options(warn=1) # print warning messages as they occur

# seed for parallelization
  # https://dofuture.futureverse.org/reference/grapes-dofuture-grapes.html 
  # https://stackoverflow.com/q/15070377
set.seed(108, kind = "L'Ecuyer-CMRG") 

### To avoid this error:
# Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) : 
# The total size of the x globals exported for future expression ('FUN()') is yyy GiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
## Set max size to 10.0 GB
options(future.globals.maxSize= 10000 * 1024^2)  # 1000 * 1024^2 = 1Gb
options(future.rng.onMisuse = "ignore")

future::plan(sequential)
print(plan())
# print(availableCores())

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
seurat_obj <- readRDS(obj_path) 

figs_dir_path <- file.path(parent_dir_path_obj, figs_dir_name)

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

#### =========================================== ####
#### Find unique markers per cluster ####
#### =========================================== ####
## Use normalized counts for DGE (Luecken & Theis 2019)
DefaultAssay(seurat_obj) <- "RNA" 

# # Use normalized counts to perform differential gene analysis
# seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE) # default: LogNormalize
# saveRDS(seurat_obj, file = file.path(subdir, paste0(filename, "_int.rds")))

print("Finding Markers...")

Idents(seurat_obj) <- cluster_col
idents <- sort(unique(seurat_obj@meta.data[[cluster_col]]))

future::plan(multisession, workers = min(floor(length(idents)/2),7)) # 3.6-re-int_niagara_gbm.sh
# future::plan(multisession, workers = floor(length(idents)/2)) # 3.6-re-int_niagara_healthy.sh
print(plan())

markers <- future_lapply(idents, function(x) {
  print(paste("Processing cluster:", x, "from resolution", cluster_col, "PID:", Sys.getpid(), "Workers Available:", nbrOfWorkers()))
  cluster_markers <- FindMarkers(
    seurat_obj, 
    ident.1 = x, 
    min.pct = 0.25, 
    logfc.threshold = 0.1)
  cluster_markers$cluster <- x
  cluster_markers$colname <- cluster_col
  cluster_markers
}, future.seed = TRUE) %>% bind_rows

future::plan(sequential)
print(plan())

write.csv(markers, file = file.path(figs_dir_path, paste0(filename, "_markers_",cluster_col,".csv")), row.names=TRUE)

#### End of Script
sessionInfo()