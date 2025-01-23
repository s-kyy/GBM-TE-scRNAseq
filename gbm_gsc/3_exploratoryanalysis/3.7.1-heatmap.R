#!usr/bin/env Rscript

# .libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
# .libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 2 filepath must be supplied: [xxx.rds] [markers.csv] [cluster_res] [fig_dir_name]", call.=FALSE)
} else if (length(args)>=2) {
  
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2])) { 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    # parent_dir_name_obj <- basename(parent_dir_path_obj)
    marker_path <- args[2]
  } else {
    stop("Filepaths provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>2) {
    cluster_resolution_num <- as.numeric(args[3])
    cluster_resolution <- paste0("integrated_snn_res.",cluster_resolution_num)
  } else {
    cluster_resolution_num <- 0.6
    cluster_resolution <- paste0("integrated_snn_res.",cluster_resolution_num)
  }
  if (length(args)>3) {
    figs_dir_name <- args[4]
  } else {
    figs_dir_name <- "figs_heatmap"
  }

} else {
    stop("Error in calling R or filepaths. Exiting...", call.=FALSE)
}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 
library(dplyr)
library(viridis)

set.seed(108)
options(warn=1) #print warning messages as they occur

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
seurat_obj <- readRDS(obj_path) 
cluster_markers <- read.csv(marker_path)
size <- 5

figs_dir_path <- file.path(parent_dir_path_obj,figs_dir_name)

ifelse(!dir.exists(file.path(figs_dir_path)),
        dir.create(file.path(figs_dir_path),recursive=T),
        "Directory Exists")

#### =========================================== ####
#### Generate Heatmap ####
#### =========================================== ####

# gather top10 genes per cluster
top10 <- cluster_markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::arrange(p_val.bonf, desc(avg_log2FC), desc(avg.exp)) %>% 
    slice_head(n = 10) %>%
    ungroup() 

# Scale and generate Heatmap
DefaultAssay(seurat_obj) <- "integrated"
Idents(seurat_obj) <- cluster_resolution

condition <- top10$gene %in% rownames(seurat_obj@assays$integrated@scale.data)

p <- DoHeatmap(seurat_obj, features = unique(top10$gene[condition])) 
ggsave(file.path(figs_dir_path, paste0("heatmap",cluster_resolution,".tiff")),
  plot = p, units="in", width=size*3.8, height=size*2.5, dpi=300, compression = 'lzw')

### End of Script
sessionInfo()