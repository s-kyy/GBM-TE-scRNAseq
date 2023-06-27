#!usr/bin/env Rscript

# R version 4.0
# Samantha Yuen
# 2023-03-07

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
library(tidyverse)

## Import data and set ident
ge_mtx <- Read10X("/home/samkyy/projects/def-ytanaka/samkyy/gete-gbm/data/2023-03-02_aggr_gbm/aggr_ge/filtered_feature_bc_matrix")
te_mtx <- Read10X("/home/samkyy/projects/def-ytanaka/samkyy/gete-gbm/data/2023-03-02_aggr_gbm/aggr_te/filtered_feature_bc_matrix")
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

ge_original <- readRDS("~/projects/def-ytanaka/common/te-gbm_proj/analysis/ge_gbmscIntUmap-subtypes.rds")
print("Imported original ge object as ge original")
gte_original <- readRDS("~/projects/def-ytanaka/common/te-gbm_proj/analysis/gte_gbmscIntUmap-subtypes.rds")
print("Imported original gte object as gte_original")

# Find intersecting cell names between genes and retrotransposon datasets
cellnames <- intersect(colnames(ge_mtx), colnames(te_mtx))
print(paste0("Number of intersecting cell names", length(cellnames)))

# Merge sparse matrices with intersecting cell names
inter <- rbind(ge_mtx[,cellnames], te_mtx[,cellnames])

# Merge with GSC sparse matrices
print(paste0("SC1 = ",unique(ge_original$sampleCombined[which(ge_original$neuro.ident == "SC1")])))
print(paste0("SC2 = ",unique(ge_original$sampleCombined[which(ge_original$neuro.ident == "SC2")])))
print(paste0("SC3 = ",unique(ge_original$sampleCombined[which(ge_original$neuro.ident == "SC3")])))
ge_original_GSC <- subset(ge_original, subset = (neuro.ident == "SC1" | neuro.ident == "SC2" | neuro.ident == "SC3"))
ge_original_GSC@meta.data <- ge_original_GSC@metadata %>%
    select(c("orig.ident", "nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio", "sample"))
rm(ge_original)
gc()

print(paste0("SC1 = ",unique(gte_original$sampleCombined[which(gte_original$neuro.ident == "SC1")])))
print(paste0("SC2 = ",unique(gte_original$sampleCombined[which(gte_original$neuro.ident == "SC2")])))
print(paste0("SC3 = ",unique(gte_original$sampleCombined[which(gte_original$neuro.ident == "SC3")])))
gte_original_GSC <- subset(gte_original, subset = (neuro.ident == "SC1" | neuro.ident == "SC2" | neuro.ident == "SC3"))
gte_original_GSC@meta.data <- gte_original_GSC@metadata %>%
    select(c("orig.ident", "nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio", "sample"))
rm(gte_original)
gc()

# Create Seurat Objects
gte <- CreateSeuratObject(counts = inter, 
                            project = "Single cell GBM Analysis -- hg38+TE",
                            names.field = 2, names.delim = "-",
                            min.features = 100,
                            min.cells = 1)
gte@meta.data <- gte@metadata %>%
    select(c("orig.ident", "nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio", "sample"))
gte_combined <- merge(gte, y=gte_original_GSC, 
                        add.cell.ids = c("GBM", "GSC"),
                        project = "Single cell GBM+GSC Analysis -- hg38+TE")
rm(te_mtx)
rm(inter)
rm(gte)
rm(gte_original_GSC)
gc()
print("Metadata of gte_combined")
print(head(gte_combined@meta.data))

ge <- CreateSeuratObject(counts = ge_mtx[,cellnames], 
                            project = "Single cell GBM Analysis -- hg38",
                            names.field = 2, names.delim = "-",
                            min.features = 100,
                            min.cells = 1)
ge_combined <- merge(ge, y=ge_original_GSC, 
                        add.cell.ids = c("GBM", "GSC"),
                        project = "Single cell GBM+GSC Analysis -- hg38")
rm(ge_mtx)
rm(ge)
rm(ge_original_GSC)
gc()
print("Metadata of ge_combined")
print(head(ge_combined@meta.data))

## Calculate quality control metrics 
# number of genes per UMI in GE+TE dataset
gte_combined$log10GenesPerUMI <- log10(gte_combined$nFeature_RNA) / log10(gte_combined$nCount_RNA)
head(gte_combined$log10GenesPerUMI)
ge_combined$log10GenesPerUMI <- log10(ge_combined$nFeature_RNA) / log10(ge_combined$nCount_RNA)
head(ge_combined$log10GenesPerUMI)

# mitochondrial ratio
gte_combined$mitoRatio <- PercentageFeatureSet(object = gte_combined, pattern = "^MT-")
gte_combined$mitoRatio <- gte_combined@meta.data$mitoRatio / 100
ge_combined$mitoRatio <- PercentageFeatureSet(object = ge_combined, pattern = "^MT-")
ge_combined$mitoRatio <- ge_combined@meta.data$mitoRatio / 100

## save objects
saveRDS(gte_combined, file="gte_combined.rds")
print(paste0("Exported initial gte seurat object"))
saveRDS(ge_combined, file="ge_combined.rds")
print(paste0("Exported initial ge seurat object"))

sessionInfo()
print("Script complete. Exiting")