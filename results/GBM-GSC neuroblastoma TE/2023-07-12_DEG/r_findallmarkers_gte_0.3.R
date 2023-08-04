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
gte <- readRDS("~/projects/def-ytanaka/common/te-gbm_proj/analysis/gte_gbmscIntUmap-subtypes.rds")
colnames(gte@meta.data)
head(gte$integrated_snn_res.0.3)
head(Idents(gte))
DefaultAssay(gte) <- "integrated"

## DEG of resolution 0.3 clusters
Idents(gte) <- "integrated_snn_res.0.3"

## run function
gte_markers <- FindAllMarkers(gte, only.pos = FALSE, min.pct = 0.05, logfc.threshold = 0.25)
    # min.diff.pct = -Inf (default to see all genes regardless of how little differences there are between groups)

## export to csv file
write.csv(gte_markers, "./gte_markers_0.3.csv", row.names=TRUE)

sessionInfo()
print("Script complete. Exiting")
