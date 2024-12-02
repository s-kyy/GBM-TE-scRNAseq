#!usr/bin/env Rscript

#### =========================================== ####
#### Verify Args ####
#### =========================================== ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
# if (length(args)<3) {
#   stop("At least 2 filepaths names must be supplied: [dataset/ge.rds] [dataset/gte.rds] [dataset/sample_counts.csv]. Optionally include path to an output folder [output_path]", call.=FALSE)
# } else if (length(args)>=3) {
if (length(args)<2) {
  stop("At least 2 filepaths names must be supplied: [dataset/ge.rds] [dataset/gte.rds]. Optionally include path to an output folder [output_path]", call.=FALSE)
} else if (length(args)>=2) {
  # verify filepaths
  # if (file.exists(args[1]) & file.exists(args[2]) & file.exists(args[3])){ 
  if (file.exists(args[1]) & file.exists(args[2])){ 
    path_ge <- args[1]  
    path_gte <- args[2] 
    # path_qc_table <- args[3]
    ge_filename <- basename(path_ext_remove(path_ge))
    gte_filename <- basename(path_ext_remove(path_gte))
    parent_dir_path <- dirname(path_ge)
    parent_dir_name <- basename(parent_dir_path)
  } else {
    stop("one or more filepaths do not exist. Closing script...", call=FALSE)
  }
  # Optional arguements
  # if (length(args) == 4 & dir.exists(args[4])) {
    # parent_dir_path <- args[4]
  if (length(args) == 3 & dir.exists(args[3])) {
    parent_dir_path <- args[3]
    parent_dir_name <- basename(parent_dir_path)
  } else if (length(args) == 3 & !dir.exists(args[3])) {
    cat("Output directory does not exist, creating new output directory...")
    dir.create(parent_dir_path, recursive=T)
  }
} 

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 

set.seed(108)

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
ge <- readRDS(path_ge) 
gte <- readRDS(path_gte) 
# qc.table <- read.csv(path_qc_table, header = TRUE)
# rownames(qc.table) <- qc.table[,1]

#### =========================================== ####
#### Filter Low Quality Cells ####
#### =========================================== ####
# Keep genes expressed in >= 3 nuclei 
gte_temp <- CreateSeuratObject(
    counts=gte$RNA@counts,
    project="reference-genes_retrotransposons", 
    min.cells = 3) # modified from createSeuratObj.R
# gte.orig_id <- gte@meta.data$orig.ident
# names(gte.orig_id) <- rownames(gte@meta.data)

ge_temp <- CreateSeuratObject(
    counts=ge$RNA@counts,
    project="reference-genes",
    min.cells = 3) # modified from createSeuratObj.R

print(paste("Original GTE object dimensions (gene x cell):", dim(gte)[1], dim(gte)[2], "down to", dim(gte_temp)[1], dim(gte_temp)[2]))
print(paste("Original GE object dimensions (gene x cell):", dim(ge)[1], dim(ge)[2], "down to", dim(ge_temp)[1], dim(ge_temp)[2]))
barcode.intersect <- intersect(colnames(gte_temp), colnames(ge_temp))
gte <- subset(gte, cells=barcode.intersect)
ge <- subset(ge, cells=barcode.intersect)
rm("gte_temp", "ge_temp")

ge_list <- SplitObject(ge, split.by="sample")
ge_list <- ge_list[order(names(ge_list))]
rm("ge")
temp_list <- vector(mode='list', length=length(ge_list))

for (i in 1:length(ge_list)) {
  print(paste("GE Sample",names(ge_list)[i],"ncells:", dim(ge_list[[i]])[2]))
  
  temp_list[[i]] <- subset(
    x=ge_list[[i]], 
    subset= (nUMI >= 500) & # remove cells with transient background level reads
      (log10GenesPerUMI > 0.80) & 
      (mitoRatio < 0.05) # nuclei should have no mitochondrial genes (with some margin)
  ) 
  print(paste("GE Sample",names(temp_list)[i],"ncells after removing background cells:", dim(temp_list[[i]])[2]))
  # # remove cells over 3 MADs in nUMI as they are likely doublets
  # # remove cells over 3 MADs in nGene as they are likely doublets
  # temp_list[[i]] <- subset(
  #   x=temp_list[[i]], 
  #   subset= (nUMI < 3 * qc.table[names(ge_list)[i], "numi.mad"]) & 
  #           (nGene < 3 * qc.table[names(ge_list)[i], "ngene.mad"]) 
  # ) 
  # print(paste("GE Sample",names(temp_list)[i],"ncells after removing doublets/triplets:", dim(temp_list[[i]])[2]))
}

filt_ge <- merge(temp_list[[1]],temp_list[-c(1)])
rm("temp_list", "ge_list")

gte_list <- SplitObject(gte, split.by="sample")
gte_list <- gte_list[order(names(gte_list))]
rm("gte")
temp_list <- vector(mode='list', length=length(gte_list))

for (i in 1:length(gte_list)) {
  print(paste("GTE Sample",names(gte_list)[i],"ncells:", dim(gte_list[[i]])[2]))
  temp_list[[i]] <- subset(
    x=gte_list[[i]], 
    subset= (nUMI >= 500) & # remove cells with transient background level reads
      (log10GenesPerUMI > 0.80) & 
      (mitoRatio < 0.05) # nuclei should have no mitochondrial genes (with some margin)
  ) 
  print(paste("GTE Sample",names(temp_list)[i],"ncells after removing background cells:", dim(temp_list[[i]])[2]))
  
  # # remove cells over 3 MADs in nUMI as they are likely doublets
  # # remove cells over 3 MADs in nGene as they are likely doublets
  # temp_list[[i]] <- subset(
  #   x=temp_list[[i]], 
  #   subset= (nUMI < 3 * qc.table[names(gte_list)[i], "numi.mad"]) &   
  #           (nGene < 3 * qc.table[names(gte_list)[i], "ngene.mad"]) 
  # ) 
  # print(paste("GTE Sample",names(temp_list)[i],"ncells after removing doublets/triplets:", dim(temp_list[[i]])[2]))
}

filt_gte <- merge(temp_list[[1]],temp_list[-c(1)])
rm("temp_list", "gte_list")

print(paste("Filtered GTE ncells:", dim(filt_gte)[2]))
print(paste("Filtered GE ncells:", dim(filt_ge)[2]))
barcode.intersect <- intersect(colnames(filt_gte), colnames(filt_ge))
filt_gte <- subset(filt_gte, cells=barcode.intersect)
filt_ge <- subset(filt_ge, cells=barcode.intersect)
print(paste("Filtered & matched GTE object ncells :", dim(filt_gte)[2]))
print(paste("Filtered & matched GE object ncells:", dim(filt_ge)[2]))

saveRDS(filt_gte, file = file.path(parent_dir_path, paste0(gte_filename, "_qc.rds")) )
saveRDS(filt_ge,  file = file.path(parent_dir_path, paste0(ge_filename, "_qc.rds")) )
cat("Filtered & saved seurat objects to",parent_dir_path,"\n")

print("Completed filtering, exiting script.")

#### End of Script ####
sessionInfo()
