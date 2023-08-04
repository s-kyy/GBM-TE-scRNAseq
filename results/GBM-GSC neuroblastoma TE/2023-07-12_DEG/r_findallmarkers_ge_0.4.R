#!usr/bin/env Rscript

# R version 4.0

set.seed(108)

## Import arguements
source("~/scratch/gete-gbm/bin/util.R")
setwd("~/scratch/temp/")
resultsPath <- "~/scratch/temp"
print(paste0("Results path:", getwd()))
mkdirToday()
print(paste0("Output path:", getwd()))

## Set library path
.libPaths(c("/scratch/samkyy/gete-gbm/renv/library/R-4.0/x86_64-pc-linux-gnu",
            "/tmp/RtmpJsRC8Z/renv-system-library",
            .libPaths()))
.libPaths()

## Import libraries
library(Seurat)
library(Matrix)

## Import data and set ident
ge <- readRDS("~/projects/def-ytanaka/common/te-gbm_proj/analysis/ge_gbmscIntUmap-subtypes.rds")
colnames(ge@meta.data)
head(ge$integrated_snn_res.0.4)
head(Idents(ge))
DefaultAssay(ge) <- "integrated"

Idents(ge) <- "integrated_snn_res.0.4"

## run function
ge_markers <- FindAllMarkers(ge, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)

## export to csv file
write.csv(ge_markers, "./ge_markers_0.4.csv", row.names=TRUE)

sessionInfo()
print("Script complete. Exiting")