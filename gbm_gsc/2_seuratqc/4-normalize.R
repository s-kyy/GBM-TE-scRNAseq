#!usr/bin/env Rscript

#### =========================================== ####
#### Verify Args ####
#### =========================================== ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("At least 1 filepaths names must be supplied: [dataset/ge.rds]. Optionally include path to an output folder [output_path]", call.=FALSE)
} else if (length(args)>=1) {
  # verify filepaths
  if (file.exists(args[1])){ 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    path_to_object <- dirname(obj_path)
    parent_dir_name <- basename(path_to_object)
  } else {
    stop("one or more filepaths do not exist. Closing script...", call=FALSE)
  }
  # Optional arguements
  if (length(args) == 2 & dir.exists(args[2])) {
    path_to_object <- args[2]
    parent_dir_name <- basename(path_to_object)
  } else if (length(args) == 2 & !dir.exists(args[2])) {
    cat("Output directory does not exist, creating new output directory...")
    dir.create(path_to_object, recursive=T)
  }
} 

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
# .libPaths(c("~/R/x86_64-pc-linux-gnu-library/tcga-gbm-R4/renv/library/R-4.0/x86_64-pc-linux-gnu", .libPaths()))
# .libPaths()

set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 

set.seed(108)

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
obj <- readRDS(obj_path) 
ndims <- 25

#### =========================================== ####
#### Normalize each dataset individually ####
#### =========================================== ####
obj_temp <- CreateAssayObject(counts = obj@assays[["RNA"]]@counts)

obj_temp <- NormalizeData(obj_temp, 
                          verbose = T, 
                          normalization.method = "LogNormalize") %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData()

obj@assays[["RNA"]]@counts <- obj_temp@counts
obj@assays[["RNA"]]@scale.data <- obj_temp@scale.data

rm("obj_temp")
gc()

obj_list <- SplitObject(object=obj, split.by="orig.ident")

rm("obj")
gc()


# Change active.ident as orig.ident for GBM samples
for (i in 1:length(obj_list)) {
    obj_list[[i]]@active.ident <- factor(obj_list[[i]]@meta.data$orig.ident)
}

# Normalize the RNA reads per tumor sample via Log Normalization method

for(i in 1:length(obj_list)){
    
    print(paste0("Processing Sample Number: ",i))
    
    DefaultAssay(obj_list[[i]]) <- "RNA"
    
    obj_list[[i]] <- NormalizeData(obj_list[[i]],verbose=TRUE)
    
    obj_list[[i]] <- FindVariableFeatures(obj_list[[i]],
                                          selection.method="vst",
                                          nfeatures=2500,
                                          verbose=TRUE)
}
cat("Variable Features computed\n")

#### =========================================== ####
#### Integrate Seurat Objects together: 30min-1h30min ####
#### =========================================== ####
obj_anchors <- FindIntegrationAnchors(obj_list,dims=1:ndims)
obj_integrated <- IntegrateData(anchorset = obj_anchors,dims=1:ndims)
cat("Integration Complete")

saveRDS(obj_integrated, file = file.path(path_to_object, paste0(filename, "_integrated.rds")) )
cat("Saved integrated seurat objects to",path_to_object,"\n")

rm("obj_list")
rm("obj_anchors")
gc()

#### =========================================== ####
#### Scale, PCA, UMAP ~5min, <200Mb ####
#### =========================================== ####
# Standard workflow for visualization and clustering 
# (ScaleData, RunPCA, RunUMAP, FindNeighbors, FindClusters)
obj_integrated <- ScaleData(object = obj_integrated, verbose = FALSE)

obj_integrated@meta.data$sample <- factor(obj_integrated@meta.data$sample, 
                                          levels = sort(unique(obj_integrated@meta.data$sample)))

# Run PCA and UMAP ~20min
obj_integrated_pca <- RunPCA(object = obj_integrated, npcs = 25, verbose = FALSE)

rm("obj_integrated")
gc()

obj_integrated_pca <- RunUMAP(object = obj_integrated_pca, 
                              reduction = "pca", 
                              dims = 1:ndims, 
                              umap.method = "uwot", 
                              metric = "cosine")

saveRDS(obj_integrated_pca, file = file.path(path_to_object, paste0(filename, "_integrated_umap.rds")) )
cat(paste("ScaleData, RunPCA, and RunUMAP complete.\n Saved Seurat objects to",path_to_object, "\n"))

#### End of Script ####
sessionInfo()
