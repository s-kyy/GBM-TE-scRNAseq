#!usr/bin/env Rscript

#### Parse Arguments ####
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<=3) {
  stop("At least 2 filepaths and 2 dataset names must be supplied: [dataset1/ge.rds] [dataset2/gte.rds] [dataset1] [dataset2]", call.=FALSE)
} else if (length(args)==4) {
  # verify filepaths
  if (file.exists(args[1]) & file.exists(args[2])){ 
    path_a <- args[1]  
    path_b <- args[2] 
  } else {
    stop("filepaths for dataset1 and/or dataset2 does not exist. Closing script...", call=FALSE)
  }
  dataset_a_name <- args[3]
  dataset_b_name <- args[4]
}

#### Import Packages ####
set.seed(108)

library(Seurat)
library(Matrix)
library(tidyverse)

set.seed(108)

#### Load Datasets ####
dataset_a <- readRDS(path_a) # USER INPUT
dataset_b <- readRDS(path_b) # USER INPUT

#### Merge Datasets into 1 Seurat Object ####

# Split datasets into a list of individual samples
# dataset_a.list <- SplitObject(object=dataset_a, split.by="orig.ident")
dataset_b.list <- SplitObject(object=dataset_b, split.by="orig.ident")
rm("dataset_b")
gc()

# Merge seurat objects
merged <- merge( x = dataset_a, 
                    y = dataset_b.list,
                    add.cell.ids = c( dataset_a_name,# c(rep(dataset_a_name,length(dataset_a.list)), 
                                     rep(dataset_b_name,length(dataset_b.list))), 
                    project = paste0("merged_", dataset_a_name, "_", dataset_b_name),
                    merge.data = FALSE) ## do not merge normalized datasets together
rm("dataset_a", "dataset_b.list")
gc()

# Export objects
subdir <- paste0(format(Sys.Date(), "%Y%m%d"), "_merged_", dataset_a_name, "_", dataset_b_name)
ifelse(!dir.exists(file.path(getwd(),subdir)),
        dir.create(file.path(getwd(),subdir),recursive=T),
        "Directory Exists")
saveRDS(merged, file = file.path(getwd(),subdir,paste0("merged_", basename(path_b))))
