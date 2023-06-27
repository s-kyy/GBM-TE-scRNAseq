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
gbmsc <- readRDS("~/projects/def-ytanaka/common/te-gbm_proj/analysis/gte_gbmscIntUmap-subtypes.rds")
# cds <- readRDS("/home/samkyy/projects/def-ytanaka/common/te-gbm_proj/analysis/2022_01_30_monocle/gbmsc_pseudo_clus6.rds"

#4. Create par object
head(gbmsc@meta.data)
Idents(gbmsc) <- gbmsc$`integrated_snn_res.0.3`
DefaultAssay(gbmsc) <- "integrated"
cds <- as.cell_data_set(gbmsc, group.by = gbmsc$integrated_snn_res.0.3)
rm(gbmsc)
gc()

#4.1 Group cells into clusters
cds <- cluster_cells(cds, random_seed = 108)

#4.2 Print cells to cluster
p1 <- plot_cells(cds, show_trajectory_graph = FALSE, group_label_size = 4)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE, group_label_size = 4)

size = 7
p1
ggsave("r_gbmsc_gte_PC20r03_monocleUMAPs_clusters_seed108.tiff", units="in", width=size, height=size, dpi=300, compression = 'lzw')
p2
ggsave("r_gbmsc_gte_PC20r03_monocleUMAPs_partitions_seed108.tiff", units="in", width=size, height=size, dpi=300, compression = 'lzw')

par <- subset(as.Seurat(cds), monocle3_partitions == 1)
# Convert back to monocle object and learn the graph
par <- as.cell_data_set(par)
par <- learn_graph(par)
# saveRDS(par, file="gbmsc_ge_pseudo_c7.rds")

#5. Select root cells
# root_cells_gte <- head (rownames(gbmsc@meta.data[which(FetchData(gbmsc,"integrated_snn_res.0.3") == 6),]))
# rm(gbmsc)
# gc()

#6. Order Cells

#6.1 FUNCTION: programmatically choose which cells come first in the trajectory
get_earliest_principal_node <- function(cds, clus="10"){
  # cell_ids <- which(colData(cds)[, "integrated_snn_res.0.3"] == clus)
  cell_ids <- which(colData(cds)[, "monocle3_clusters"] == clus)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

par_auto <- order_cells(par, root_pr_nodes=get_earliest_principal_node(par))

saveRDS(par_auto,paste0(resultsPath, "/gte_monocle3.rds"))
# save_monocle_objects(cds=par_auto, 
#                     directory_path=getwd(), 
#                     comment='GBMGSC analysis: hg38+TE sorted cells, root=cluster6, Stored 2023-02-13.')

p2 <- par_auto %>% plot_cells(color_cells_by = "pseudotime", 
           label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = FALSE,
           label_roots = TRUE,
           trajectory_graph_color = "black",
           graph_label_size = 3,
           cell_size = 0.5)
p2

size=4
ggsave("r_gbmsc_gte_PC20r03_monocleUMAPs_path.tiff", units="in", width=size, height=size, dpi=300, compression = 'lzw')

sessionInfo()