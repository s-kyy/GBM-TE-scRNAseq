#!usr/bin/env Rscript

# R version 4.0
# Samantha Yuen
# 2023-03-20

set.seed(108)

## Import arguements
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("-input_ge_dir -input_te_dir -sampleIDs.csv are mandatory. Optional: -outputfilename -outputdir", call.=FALSE)
} else if (length(args)==1) {
    stop("-input_ge_dir -input_te_dir -sampleIDs.csv are mandatory. Optional: -outputfilename -outputdir", call.=FALSE)
} else if (length(args)==2) {
    stop("-input_ge_dir -input_te_dir -sampleIDs.csv are mandatory. Optional: -outputfilename -outputdir", call.=FALSE)
} else if (length(args)==3) {
    args[4] <- "out.rds"
    args[5] <- getwd()
} else if (length(args)==4) {
    args[5] <- getwd()
}
if (!dir.exists(args[1])) {stop("input_ge directory does not exist", call.=FALSE)}
if (!dir.exists(args[2])) {stop("input_te directory does not exist", call.=FALSE)}
if (!file.exists(args[3])) {stop("sampleIDs csv does not exist", call.=FALSE)}
print("Arguements passed:")
print(args)

source("~/scratch/gete-gbm/bin/util.R")
resultsPath <- args[5]
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
library(tidyverse)

## Import data and set ident
# EG: /home/samkyy/projects/def-ytanaka/samkyy/gete-gbm/data/2023-03-06_aggr_normal/aggr_ge/outs/count/filtered_feature_bc_matrix
ge_mtx <- Read10X(args[1])
te_mtx <- Read10X(args[2])
print("Dimensions of hg38 annotated matrix: ")
print(dim(ge_mtx))
print("First ten cell names ")
print(colnames(ge_mtx)[1:10])
print("First ten genes")
print(rownames(ge_mtx)[1:10])
print("Dimensions of TE annotated matrix: ")
print(dim(te_mtx))
print("First ten cell names")
print(colnames(te_mtx)[1:10])
print("First ten genes")
print(rownames(te_mtx)[1:10])

samples <- as.list(read.csv(args[3], header = FALSE))
samples <- samples$V1
print(samples)
# Find intersecting cell names between genes and retrotransposon datasets
cellnames <- intersect(colnames(ge_mtx), colnames(te_mtx))
print(paste0("Number of intersecting cell names", length(cellnames)))

# Merge sparse matrices with intersecting cell names
inter <- rbind(ge_mtx[,cellnames], te_mtx[,cellnames])

# Create Seurat Objects
gte <- CreateSeuratObject(counts = inter, 
                            project = "Single cell GBM Analysis -- hg38+TE",
                            names.field = 2, names.delim = "-",
                            min.features = 100,
                            min.cells = 1)
rm(te_mtx)
rm(inter)
gc()

ge <- CreateSeuratObject(counts = ge_mtx[,cellnames], 
                            project = "Single cell GBM Analysis -- hg38",
                            names.field = 2, names.delim = "-",
                            min.features = 100,
                            min.cells = 1)
rm(ge_mtx)
gc()


## Calculate quality control metrics 
# number of genes per UMI in GE+TE dataset
gte$log10GenesPerUMI <- log10(gte$nFeature_RNA) / log10(gte$nCount_RNA)
head(gte$log10GenesPerUMI)
ge$log10GenesPerUMI <- log10(ge$nFeature_RNA) / log10(ge$nCount_RNA)
head(ge$log10GenesPerUMI)

# mitochondrial ratio
gte$mitoRatio <- PercentageFeatureSet(object = gte, pattern = "^MT-")
gte$mitoRatio <- gte@meta.data$mitoRatio / 100
ge$mitoRatio <- PercentageFeatureSet(object = ge, pattern = "^MT-")
ge$mitoRatio <- ge@meta.data$mitoRatio / 100
print("\n")
print("Metadata of combined seurat objects")
gte@meta.data %>% glimpse()
print("\n")
ge@meta.data %>% glimpse()
print("\n")

# Clean metadata
ge_meta <- ge@meta.data
gte_meta <- gte@meta.data

ge_meta$sampleID = ""
gte_meta$sampleID = ""

ge_meta %>% glimpse()
gte_meta %>% glimpse()

for (i in 1:length(samples)) {
    ge_meta$sampleID[which(ge_meta$orig.ident == i)] <- samples[i]
    gte_meta$sampleID[which(gte_meta$orig.ident == i)] <- samples[i]
}

print(unique(ge_meta$sampleID))
print(unique(gte_meta$sampleID))

ge@meta.data <- ge_meta
gte@meta.data <- gte_meta

## save objects
saveRDS(gte, file=paste0("gte_", args[4]))
print(paste0("Exported initial gte seurat object"))
saveRDS(ge, file=paste0("ge_", args[4]))
print(paste0("Exported initial ge seurat object"))

sessionInfo()
print("Script complete. Exiting")