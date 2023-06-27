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
ggsave("r_gbmsc_gte_PC20r03_monocleUMAPs_clusters_seed108.tiff", units="in", width=size, height=size, dpi=300, compression = 'lzw')
p2
ggsave("r_gbmsc_gte_PC20r03_monocleUMAPs_partitions_seed108.tiff", units="in", width=size, height=size, dpi=300, compression = 'lzw')

sessionInfo()