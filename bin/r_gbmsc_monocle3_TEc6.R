#!usr/bin/env Rscript

## Import arguements
source("~/scratch/gete-gbm/bin/util.R")
resultsPath <- "~/scratch/gete-gbm/results/"

# args <- commandArgs(trailingOnly = TRUE)
# if ((length(args) == 0)) {

# 	stop("3 parameter values are required: [path to seurat object] [filepath to output directory]  \nSeurat object must be a .rds file", call= FALSE)

# } else if (length(args) == 1) {
    
#     if (file.exists(args[1])) { dataPath <- args[1] 
#     } else { stop("Path to data does not exist, please try again ", call = FALSE)}
#     ## resultsPath is "~/scratch/temp/" by default

# } else {
    
#     if (file.exists(args[1])) { dataPath <- args[1]
#     } else { stop("Path to data does not exist, please try again ", call = FALSE)}

#     if (dir.exists(dirname(args[2]))) { resultsPath <- args[2]
#     } else { stop("Output directory for .rds file does not exist, please try again ", 
#             call = FALSE)}
# }

# cat("Variables Imported\n")

## Set library path
.libPaths(c("/scratch/samkyy/gete-gbm/renv/library/R-4.0/x86_64-pc-linux-gnu",
    "/tmp/RtmpJsRC8Z/renv-system-library",
    .libPaths()))
.libPaths()
cat("\n")

## Load libraries
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

# par1 <- readRDS(dataPath)
# gbmsc <- readRDS("~/scratch/gete-gbm/results/2021-09-02/gte_gbmscIntUmap-subtypes.rds")
par1 <- readRDS("~/scratch/gete-gbm/results/2022-01-30/gbmsc_pseudo_clus6.rds")

p1 <- par1 %>% plot_cells(
        label_groups_by_cluster = FALSE,
        label_leaves = TRUE,
        label_branch_points = FALSE,
        graph_label_size = 3,
        # trajectory_graph_color = "white"
        )
p1

## Select root cells
    ##GTE_PC20_r0.3: Cluster 14

root_cells_gte <- head(rownames(gbmsc@meta.data[which(FetchData(gbmsc,"integrated_snn_res.0.3") == 6),]))

## Order Cells
par1_c6 <- order_cells(par1, root_cells = root_cells_gte)

p2 <- par1_c6 %>% plot_cells(color_cells_by = "pseudotime", 
           label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "black",
           graph_label_size = 5,
           cell_size = 0.5)
p2

## FUNCTION: programmatically choose which cells come first in the trajectory
get_earliest_principal_node <- function(cds, clus="6"){
  cell_ids <- which(colData(cds)[, "integrated_snn_res.0.3"] == clus)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

par1_auto <- order_cells(par1, root_pr_nodes=get_earliest_principal_node(par1))

p2 <- par1_auto %>% plot_cells(color_cells_by = "pseudotime", 
           label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = FALSE,
           label_roots = TRUE,
           trajectory_graph_color = "black",
           graph_label_size = 3,
           cell_size = 0.5)
p2

setwd(resultsPath)
getwd()
mkdirToday()
# [1] "Current working directory: /scratch/samkyy/gete-gbm/results"
# [1] "New working directory: /scratch/samkyy/gete-gbm/results/2022-02-09"

ggsave("r_gbmsc_PC20r03_monocleUMAPs_path.tiff", units="in", width=size, height=size, dpi=300, compression = 'lzw')