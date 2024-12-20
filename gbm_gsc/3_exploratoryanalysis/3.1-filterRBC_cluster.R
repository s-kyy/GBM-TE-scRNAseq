#!usr/bin/env Rscript

#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 1 filepath must be supplied: [xxx.rds] [sample_name]", call.=FALSE)
} else if (length(args)<2) {
  # verify filepaths
  if (file.exists(args[1])){ 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    sample_name <- ""
  } else {
    stop("Filepath provided does not exist. Exiting...", call.=FALSE)
  }
} else if (length(args)<=2) {
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    sample_name <- args[2]
}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
.libPaths(c("~/scratch/tcga-gbm-R4.1-lib/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()

set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 

set.seed(108)

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
seurat.obj <- readRDS(obj_path) 
size    <- 5

subdir <- file.path(getwd(), paste0(format(Sys.Date(), "%Y%m%d"), "_", sample_name,"_", filename))

ifelse(!dir.exists(file.path(subdir)),
        dir.create(file.path(subdir),recursive=T),
        "Directory Exists")
figs_dir_path <- file.path(subdir, "figs_filt")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

## ========================================= ##
## Filter RBC cells ##
## ========================================= ##

filtered_obj <- subset(
  x=seurat.obj, 
  subset= 
  (HBA1 < 1) &
  (HBA2 < 1) &
  (HBB < 1) &
  (HBG2 < 1) &
  (HBZ < 1) & 
  (HBM < 1) & 
  (HBD < 1) & 
  (HBE1 < 1) &
  (HBG1 < 1) &
  (HBQ1 < 1)  
) 
print(paste(dim(seurat.obj)[2] - dim(filtered_obj)[2], "cells were removed from object"))

# Find variable features in the subset using the RNA assay
# Run ScaleData on the integrate assay on the new set of variable features
# Run PCA on the integrated assay using the new set of variable features
# Run FindNeighbors and FindClusters using the new PC dimensions
# source: https://github.com/satijalab/seurat/issues/5532

#### =========================================== ####
#### Normalize & Find Variable Features each dataset individually ####
#### =========================================== ####

# Correct for inter-sample variability 
# (i.e. mainly from different tumour donors) 
print(paste("unique samples to split by:", unique(seurat.obj$sample)))
obj_list <- SplitObject(object=seurat.obj, split.by="sample")
print(paste("Length of list", length(obj_list)))
rm("obj")
gc()

# Normalize the RNA reads per tumor sample via Log Normalization method
p_list <- vector(mode="list",length(obj_list))
names(p_list) <- sort(names(obj_list))

for(i in 1:length(obj_list)){
    print(paste0("Processing Sample: ", obj_list[[i]]$sample[1]))
    DefaultAssay(obj_list[[i]]) <- "RNA"
    # default: LogNormalize, save to : $RNA@data
    obj_list[[i]] <- NormalizeData(obj_list[[i]], verbose = TRUE) 
    obj_list[[i]] <- FindVariableFeatures(obj_list[[i]], 
                              selection.method="vst",
                              nfeatures=2500, #default=2000
                              verbose = TRUE)
    #Create Variable Feature Plots
    top10 <- head(VariableFeatures(obj_list[[i]]), 10)
    p <- VariableFeaturePlot(obj_list[[i]])
    p <- LabelPoints(plot = p, points = top10, repel = TRUE)
    ggsave(file.path(figs_dir_path, paste0(filename,obj_list[[i]]$sample[1],"_varfeatgenes.tiff")), 
          plot = p, units="in", width=size*1.1, height=size*0.8, dpi=300, compression = 'lzw')
    p_list[[obj_list[[i]]$sample[1]]] <- p
}

saveRDS(p_list, file=file.path(figs_dir_path,paste0(filename,"_allvarfeatplots.rds")))

print("Variable Features computed")

#### =========================================== ####
#### Integrate Seurat Objects together: 30min-1h30min ####
#### =========================================== ####
## Anchors across all samples will be used to integrated each sample
## This step helps reduce variation introduced across donor samples 
## (i.e. variation in tumour extraction, experimental processing of samples, 
## gender-based differences and so on)
## Stuck with CCA method to find anchors, because the datasets are a manageable size (<50K cells)
obj_anchors <- FindIntegrationAnchors(
                  object.list = obj_list, 
                  dims=1:ndims,
                  normalization.method = "LogNormalize", 
                  anchor.features = 2500, #match number in FindVariableFeatures
                  )
obj_integrated <- IntegrateData(anchorset = obj_anchors, dims=1:ndims)
saveRDS(obj_integrated, file = file.path(parent_dir_path,output_dir_name, paste0(filename, "_integrated.rds")) )
print("Integration Complete")

## ========================================= ##
## Scale and Reduce Dimensionality ##
## ========================================= ##
DefaultAssay(filtered_obj) <- "integrated"
# Scale Data ~5min, <200Mb
filtered_obj <- FindVariableFeatures(filtered_obj, 
                              selection.method="vst",
                              nfeatures=2500, #default=2000
                              verbose = TRUE)
filtered_obj <- ScaleData(object = filtered_obj, verbose = FALSE) 

filtered_obj@meta.data$sample <- factor(filtered_obj@meta.data$sample, 
                                          levels = sort(unique(filtered_obj@meta.data$sample)))
# Run PCA and UMAP ~20min
filtered_obj <- RunPCA(filtered_obj, 
                          npcs = ndims, 
                          verbose = FALSE)

filtered_obj <- RunUMAP(filtered_obj, 
                          reduction = "pca", 
                          dims = 1:ndims, 
                          umap.method = "uwot", 
                          metric = "cosine")
saveRDS(filtered_obj, file = file.path(subdir, paste0(filename, "_filt_umap.rds")) )
print(paste("ScaleData, RunPCA, and RunUMAP completed in:", subdir))

#### End of Script ####
sessionInfo()
