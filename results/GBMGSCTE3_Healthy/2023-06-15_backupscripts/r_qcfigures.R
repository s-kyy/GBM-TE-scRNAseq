#!usr/bin/env Rscript

# R version 4.0
# Samantha Yuen
# 2023-03-08

set.seed(108)

## Import arguements
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("-seurat-objects are mandatory. Optional: -outputfilename -outputdir", call.=FALSE)
} else if (length(args)==1) {
    args[2] <- "out"
    args[3] <- getwd()
} else if (length(args)==2) {
    args[3] <- getwd()
}
if (!file.exists(args[1])) {stop("ge_object directory does not exist", call.=FALSE)}
if (!dir.exists(args[4])) {stop("outputdir does not exist", call.=FALSE)}
print("Arguements passed:")
print(args)

source("~/scratch/gete-gbm/bin/util.R")
resultsPath <- args[3]
print(paste0("Results path:", resultsPath))
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
library(ggplot2)
library(tidyverse)

## Import data
obj <- readRDS(args[1])
size = 7 # default figure size in pixels

obj@meta.data %>% glimpse()

### Cell Counts per Sample
print("Exporting Fig: Cell Counts per sample (ID)")
if ('sampleID' %in% rownames(obj@meta.data)) {
    p1 <- obj@meta.data %>%
    ggplot(aes(x=orig.ident, fill=orig.ident)) + 
    labs(fill="samples") + labs(x = "Samples", y = "Cell Count") +
    scale_x_discrete(labels = sampleID) +
    scale_fill_discrete(labels = sampleID) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    tiff(paste0("fig_",args[2],"_CellCountPerSample.tiff"),units="in", width=size*1.5, height=size, dpi=300, compression = 'lzw')
    p1
    dev.off()
}

print("Exporting Fig: Cell Counts per sample (full ID)")
if ('sample_full' %in% rownames(obj@meta.data)) {
    p1 <- obj@meta.data %>%
    ggplot(aes(x=orig.ident, fill=orig.ident)) + 
    labs(fill="samples") + labs(x = "Samples", y = "Cell Count") +
    scale_x_discrete(labels = sampleID) +
    scale_fill_discrete(labels = sampleID) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    tiff(paste0("fig_",args[2],"_CellCountPerSample(full).tiff"),units="in", width=size*1.5, height=size, dpi=300, compression = 'lzw')
    p1
    dev.off()
}

## Figures of Quality control metrics -- by sample

## Filter out cells based on quality control metrics

## save filtered objects metadata