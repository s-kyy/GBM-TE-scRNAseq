#!usr/bin/env Rscript

# R version 4.0
# Samantha Yuen
# 2023-03-17

set.seed(108)

## Import arguements
source("~/scratch/gete-gbm/bin/util.R")
setwd("~/scratch/temp/analysis")
resultsPath <- "~/scratch/temp/analysis"
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

gte <- 
ge <- 