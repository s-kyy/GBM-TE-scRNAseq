#!usr/bin/env Rscript

#### =========================================== ####
#### Verify Args ####
#### =========================================== ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
if (length(args)<1) {
  stop("At least 1 filepath names must be supplied: [dataset/ge.rds]. Optionally include path to an output folder name [output_dir_name]", call.=FALSE)
} else if (length(args)>=1) {
  # verify filepaths
  if (file.exists(args[1])){ 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path <- dirname(obj_path)
    parent_dir_name <- basename(parent_dir_path)
  } else {
    stop("filepath does not exist. Closing script...", call=FALSE)
  }
  # Optional arguements
  if (length(args) == 2) {
    output_dir_name <- args[2]
} 

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/tcga-gbm-R4/renv/library/R-4.0/x86_64-pc-linux-gnu", .libPaths()))
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
obj <- readRDS(obj_path) 
ndims <- 25
size <- 5

ifelse(!dir.exists(file.path(parent_dir_path,output_dir_name)),
        dir.create(file.path(parent_dir_path,output_dir_name),recursive=T),
        "Directory Exists")

figs_dir_path <- file.path(parent_dir_path, output_dir_name, "figs")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

#### =========================================== ####
#### Normalize & Find Variable Features each dataset individually ####
#### =========================================== ####

# Correct for inter-sample variability 
# (i.e. mainly from different tumour donors) 
print(paste("unique samples to split by:", unique(obj$sample)))
obj_list <- SplitObject(object=obj, split.by="sample")
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
## Scale and Reduce Dimensionality without regression ##
## ========================================= ##
DefaultAssay(obj_integrated) <- "integrated"
# Scale Data ~5min, <200Mb
obj_integrated <- ScaleData(object = obj_integrated, verbose = FALSE) 

obj_integrated@meta.data$sample <- factor(obj_integrated@meta.data$sample, 
                                          levels = sort(unique(obj_integrated@meta.data$sample)))
# Run PCA and UMAP ~20min
obj_integrated <- RunPCA(obj_integrated, 
                          npcs = ndims, 
                          verbose = FALSE)

obj_integrated <- RunUMAP(obj_integrated, 
                          reduction = "pca", 
                          dims = 1:ndims, 
                          umap.method = "uwot", 
                          metric = "cosine")
saveRDS(obj_integrated, file = file.path(parent_dir_path, output_dir_name, paste0(filename, "_integrated_umap.rds")) )
print(paste("ScaleData, RunPCA, and RunUMAP completed in:", parent_dir_path, output_dir_name))

#### End of Script ####
sessionInfo()
