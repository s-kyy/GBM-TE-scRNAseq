#!usr/bin/env Rscript

#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 1 filepath must be supplied: [xxx.rds]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1])){ 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    # path_to_object <- dirname(obj_path)
    # parent_dir_name <- basename(path_to_object)
  } else {
    stop("Filepath provided does not exist. Exiting...", call.=FALSE)
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
seurat.obj <- readRDS(obj_path) 
size    <- 5

subdir <- file.path(getwd(), paste0(format(Sys.Date(), "%Y%m%d"), "_", filename))

ifelse(!dir.exists(file.path(subdir)),
        dir.create(file.path(subdir),recursive=T),
        "Directory Exists")

ifelse(!dir.exists(file.path(subdir, "figs")),
        dir.create(file.path(subdir, "figs"),recursive=T),
        "Directory Exists")


#### =========================================== ####
#### Test PCA levels ####
#### =========================================== ####

# Determine percent of variation associated with each PC
pct_var_per_pc <- seurat.obj[["pca"]]@stdev / sum(seurat.obj[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cum_pct_per_pc <- cumsum(pct_var_per_pc)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
min_pc <- which(cum_pct_per_pc > 90 & pct_var_per_pc < 5)[1]
print(paste("Minimum PC that retains more than 90% variation and less than 5% variation compared to the next PC:", min_pc))

plot_df <- data.frame(dimensions = 1:length(pct_var_per_pc),
           stdev = seurat.obj[["pca"]]@stdev,
           pct_var_per_pc = pct_var_per_pc,
           cum_pct_per_pc = cum_pct_per_pc)
print(plot_df)
write.csv(plot_df, file.path(subdir, "figs", paste0(filename, "_pca.csv") ))

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
ggsave(file.path(subdir, "figs", paste0(filename,"_elbow.tiff")), 
  plot = p, units="in", width=size*0.7, height=size*0.7, dpi=300, compression = 'lzw')
print("Exported ElbowPlot")

#### =========================================== ####
#### Output UMAP plots by sample ####
#### =========================================== ####
size <- 5

p <- DimPlot(seurat.obj, reduction = "umap", group.by = "sample") 
ggsave(file.path(subdir, "figs", paste0(filename,"_UMAP-sample.tiff")), 
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')
print("Exported UMAP")

#### =========================================== ####
# Cluster: K-nearest neighbor graph
#### =========================================== ####
seurat.obj <- FindNeighbors(seurat.obj,dims=1:min_pc,reduction="pca")
seurat.obj <- FindClusters(seurat.obj, resolution = 0.3)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.4)
saveRDS(seurat.obj, file = file.path(subdir, paste0(filename, "_clustered.rds")))
cat("KNN clustering complete. Saved seurat objects to", subdir,"\n")

# Add FindAllMarkers function
DefaultAssay(seurat.obj) <- "integrated"
Idents(seurat.obj) <- "integrated_snn_res.0.3"

markers <- FindAllMarkers(
  seurat.obj, 
  only.pos = FALSE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25)

write.csv(markers, file = file.path(subdir, paste0(filename, "_markers_0.3.csv")), row.names=TRUE)
rm("markers")

## DEG of resolution 0.4 clusters
Idents(seurat.obj) <- "integrated_snn_res.0.4"

## run function
markers <- FindAllMarkers(
  seurat.obj, 
  only.pos = FALSE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25)

write.csv(markers, file = file.path(subdir, paste0(filename, "_markers_0.4.csv")), row.names=TRUE)
rm("markers")

#### =========================================== ####
#### Output UMAP plots by cluster ####
#### =========================================== ####
p <- DimPlot(seurat.obj, reduction = "umap", group.by = "cluster") 
ggsave(file.path(subdir, "figs", paste0(filename,"_UMAP-cluster.tiff")),
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

cat("Seurat object processing complete. Saved to ",subdir,"\n")

#### End of Script ####
sessionInfo()
