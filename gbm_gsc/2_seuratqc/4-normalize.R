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
    parent_dir_path <- dirname(obj_path)
    parent_dir_name <- basename(parent_dir_path)
  } else {
    stop("one or more filepaths do not exist. Closing script...", call=FALSE)
  }
  # Optional arguements
  if (length(args) == 2 & dir.exists(args[2])) {
    parent_dir_path <- args[2]
    parent_dir_name <- basename(parent_dir_path)
  } else if (length(args) == 2 & !dir.exists(args[2])) {
    print("Output directory does not exist, creating new output directory...")
    dir.create(parent_dir_path, recursive=T)
  }
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

size    <- 5

figs_dir_path <- file.path(parent_dir_path, "figs")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

#### =========================================== ####
#### Normalize each dataset individually ####
#### =========================================== ####

# Temporarily create new seurat object required when 
# seurat objects were merged (quick solution to a bug). 
obj_temp <- CreateAssayObject(counts = obj@assays[["RNA"]]@counts)

obj_temp <- NormalizeData(obj_temp, 
                          verbose = T, 
                          normalization.method = "LogNormalize") %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData()

obj@assays[["RNA"]]@counts <- obj_temp@counts
obj@assays[["RNA"]]@scale.data <- obj_temp@scale.data # not used for preprocessing

rm("obj_temp")
gc()

# Correct for inter-sample variability 
# (i.e. mainly from different tumour donors) 
obj_list <- SplitObject(object=obj, split.by="orig.ident")
rm("obj")
gc()

p_list <- vector(mode="list",length(obj_list))
names(p_list) <- sort(names(obj_list))

# Normalize the RNA reads per tumor sample via Log Normalization method
obj_list <- lapply(X = obj_list, FUN = function(x) {
    print(paste0("Processing Sample: ", x$sample[1]))
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, verbose = TRUE) # default: LogNormalize
    x <- FindVariableFeatures(x, 
                              selection.method="vst"
                              nfeatures=2500, #default=2000
                              verbose = TRUE)
    #Create Variable Feature Plots
    top10 <- head(VariableFeatures(x), 10)
    p <- VariableFeaturePlot(x)
    p <- LabelPoints(plot = p, points = top10, repel = TRUE)
    ggsave(file.path(figs_dir_path, paste0(filename,x$sample[1],"_varfeatgenes.tiff")), 
          plot = p, units="in", width=size*1.1, height=size*0.8, dpi=300, compression = 'lzw')
    y[[x$sample[1]]] <- p
})
saveRDS(p_list, file=file.path(figs_dir_path,paste0(filename,"_allvarfeatplots.rds")))

print("Variable Features computed")

#### =========================================== ####
#### Integrate Seurat Objects together: 30min-1h30min ####
#### =========================================== ####
## Anchors across all samples will be used to integrated each sample
## This step helps reduce variation introduced across donor samples 
## (i.e. variation in tumour extraction, experimental processing of samples, 
## gender-based differences and so on)
obj_anchors <- FindIntegrationAnchors(obj_list,dims=1:ndims)
obj_integrated <- IntegrateData(anchorset = obj_anchors,dims=1:ndims)
print("Integration Complete")

saveRDS(obj_integrated, file = file.path(parent_dir_path, paste0(filename, "_integrated.rds")) )
print(paste("Saved integrated seurat objects to",parent_dir_path))

rm("obj_list")
rm("obj_anchors")
gc()

#### =========================================== ####
#### Scale, PCA, UMAP ~5min, <200Mb ####
#### =========================================== ####
obj_integrated <- ScaleData(object = obj_integrated, verbose = FALSE)
  ### add vars.to.regress = c("mitoRatio")!!!
  ### https://github.com/satijalab/seurat/discussions/4259

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

saveRDS(obj_integrated_pca, file = file.path(parent_dir_path, paste0(filename, "_integrated_umap.rds")) )
print(paste("ScaleData, RunPCA, and RunUMAP complete.\n Saved Seurat objects to",parent_dir_path))

#### End of Script ####
sessionInfo()
