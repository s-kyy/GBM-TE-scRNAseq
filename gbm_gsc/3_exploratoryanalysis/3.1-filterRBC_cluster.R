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
saveRDS(filtered_obj, file = file.path(subdir, paste0(filename, "_integrated_umap.rds")) )
print(paste("ScaleData, RunPCA, and RunUMAP completed in:", subdir))

#### =========================================== ####
#### Test PCA levels ####
#### =========================================== ####
# Determine percent of variation associated with each PC
DefaultAssay(filtered_obj) <- "integrated"
pct_var_per_pc <- filtered_obj[["pca"]]@stdev / sum(filtered_obj[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cum_pct_per_pc <- cumsum(pct_var_per_pc)

# Determine which PC exhibits a cumulative percentage of variation 
# greater than 90% and variation associated with the PC is less than 5%
min_pc <- which(cum_pct_per_pc > 90 & pct_var_per_pc < 5)[1]
print(paste("Minimum PC that retains more than 90% variation and less than 5% variation compared to the next PC:", min_pc))

plot_df <- data.frame(dimensions = 1:length(pct_var_per_pc),
           stdev = seurat.obj[["pca"]]@stdev,
           pct_var_per_pc = pct_var_per_pc,
           cum_pct_per_pc = cum_pct_per_pc)
print(plot_df)
write.csv(plot_df, file.path(figs_dir_path, paste0(filename, "_pca.csv") ))

# Plot % variation to Elbowplot (modified from Seurat Elbow Plot Function)
p <- plot_df %>% ggplot(aes(x = dimensions, y = stdev)) +
    geom_point() +
    labs(x = "", y = "Standard Deviation") +
    geom_text(
      label=format(round(cum_pct_per_pc, 1), nsmall = 1), 
      nudge_x = 0.5, nudge_y = 0.5, 
      check_overlap = T,
      size=2) +
    theme_classic() 
ggsave(file.path(figs_dir_path, paste0(filename,"_elbow.tiff")), 
  plot = p, units="in", width=size*0.7, height=size*0.7, dpi=300, compression = 'lzw')
print("Exported ElbowPlot")

#### =========================================== ####
#### Output UMAP plots by sample ####
#### =========================================== ####
size <- 5

p <- DimPlot(filtered_obj, reduction = "umap", group.by = "sample") 
ggsave(file.path(figs_dir_path, paste0(filename,"_UMAP-sample.tiff")), 
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')
print("Exported UMAP by sample") 

#### =========================================== ####
#### Cluster: K-nearest neighbor graph
#### =========================================== ####
filtered_obj <- FindNeighbors(filtered_obj,dims=1:min_pc,reduction="pca")
filtered_obj <- FindClusters(filtered_obj, resolution = 0.3)
filtered_obj <- FindClusters(filtered_obj, resolution = 0.4)
saveRDS(filtered_obj, file = file.path(subdir, paste0(filename, "_clustered.rds")))
cat("KNN clustering complete. Saved seurat objects to", subdir,"\n")

#### =========================================== ####
#### Find unique markers per cluster ####
#### =========================================== ####
## Use normalized counts for DGE (Luecken & Theis 2019)
DefaultAssay(filtered_obj) <- "RNA" 
Idents(filtered_obj) <- "integrated_snn_res.0.3"

# Use normalized counts to perform differential gene analysis
filtered_obj <- NormalizeData(filtered_obj, verbose = TRUE) # default: LogNormalize

markers <- FindAllMarkers(
  filtered_obj, 
  only.pos = FALSE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25)

write.csv(markers, file = file.path(subdir, paste0(filename, "_markers_0.3.csv")), row.names=TRUE)
rm("markers")

## DEG of resolution 0.4 clusters
Idents(filtered_obj) <- "integrated_snn_res.0.4"

## run function
markers <- FindAllMarkers(
  filtered_obj, 
  only.pos = FALSE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25)

write.csv(markers, file = file.path(subdir, paste0(filename, "_markers_0.4.csv")), row.names=TRUE)
rm("markers")

#### End of Script ####
sessionInfo()
