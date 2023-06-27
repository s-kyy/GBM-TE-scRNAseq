#!/bin/bash/Rscript

# Load Librairies

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

#Call the commands below with the following structure whenever I want to auto-select a funtion from a specific package:
library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
# conflict_prefer("rowRanges", "matrixStats")
conflict_prefer("rowRanges", "MatrixGenerics") ## required in new_cell_data_set()

udunits_dir = "/home/samkyy/bin/udunits"
dyn.load(paste0(udunits_dir, "/local/lib/libudunits2.so.0"))