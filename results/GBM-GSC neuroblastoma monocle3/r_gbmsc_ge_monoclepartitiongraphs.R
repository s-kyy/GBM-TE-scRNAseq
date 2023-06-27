#!usr/bin/env Rscript

# R version 4.0

#1. Set working directory
source("~/scratch/gete-gbm/bin/util.R")
resultsPath <- "~/scratch/temp/"
setwd(resultsPath)
print(paste0("Results path:", getwd()))
mkdirToday()
print(paste0("Output path:", getwd()))

set.seed(108)


#2. Set library path
.libPaths(c("/scratch/samkyy/gete-gbm/renv/library/R-4.0/x86_64-pc-linux-gnu",
    "/tmp/RtmpJsRC8Z/renv-system-library",
    .libPaths()))
.libPaths()
cat("\n")

#2.1 Load libraries
# source("~/scratch/gete-gbm/bin/util_monocle3.R")
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
library(patchwork)
library(magrittr)
library(tidyverse)
library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("rowRanges", "MatrixGenerics") # required in new_cell_data_set()

#3. Import Dataset 
gbmsc <- readRDS("~/projects/def-ytanaka/common/te-gbm_proj/analysis/ge_gbmscIntUmap-subtypes.rds")
# par1 <- readRDS("/home/samkyy/projects/def-ytanaka/common/te-gbm_proj/analysis/2022_01_30_monocle/gbmsc_pseudo_clus6.rds")

#4. Create par object
head(gbmsc@meta.data)
Idents(gbmsc) <- gbmsc$`integrated_snn_res.0.3`
DefaultAssay(gbmsc) <- "integrated"
cds <- as.cell_data_set(gbmsc, group.by = gbmsc$integrated_snn_res.0.3)

#4.1 Group cells into clusters
cds <- cluster_cells(cds, random_seed = 108)

#4.2 Print cells to cluster
p1 <- plot_cells(cds, show_trajectory_graph = FALSE, group_label_size = 4)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE, group_label_size = 4)

size = 7
p1
ggsave("r_gbmsc_ge_PC20r03_monocleUMAPs_clusters_seed108.tiff", units="in", width=size, height=size, dpi=300, compression = 'lzw')
p2
ggsave("r_gbmsc_ge_PC20r03_monocleUMAPs_partitions_seed108.tiff", units="in", width=size, height=size, dpi=300, compression = 'lzw')

# pdf("r_gbmsc_ge_PC20r03_monocleUMAPs.pdf", width = size*2, height = size)
# grid.arrange(p1, p2, nrow = 1, top="Trajectory Analysis: Monocle Clusters vs Partitions")
# dev.off()

## Extract sparse matrix, metadata and gene labels to create monocle object
# data <- as(as.matrix(GetAssayData(gbmsc, assay = "integrated", slot = "scale.data")), 'sparseMatrix')
# meta <- data.frame(gbmsc@meta.data) %>% 
#             select(orig.ident, nCount_RNA, nFeature_RNA, sampleCombined, seurat_clusters)
# genes <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

# monocle <- new_cell_data_set(expression_data = data, cell_metadata = meta, gene_metadata = genes) 

# # normalize
# monocle <- preprocess_cds(monocle, num_dim = 100, norm_method = "size_only", pseudo_count = 0)  
# #skip batch effect removal since that is already complete with Seurat 

# # reduce dimensions 
# monocle <- reduce_dimension(monocle)


# Select Island to study the trajectory on:
# var = readline("Please enter the partition number to subset: ") 
# var = as.integer(var)
# par <- subset(as.Seurat(cds), monocle3_partitions == var)
# # Convert back to monocle object and learn the graph
# par <- as.cell_data_set(par)
# par <- learn_graph(par)
# saveRDS(par, file="gbmsc_pseudo_clus6.rds")

sessionInfo()