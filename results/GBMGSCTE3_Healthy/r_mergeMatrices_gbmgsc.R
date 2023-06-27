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
# EG: /home/samkyy/projects/def-ytanaka/samkyy/gete-gbm/data/2023-03-02_aggr_gbm/aggr_ge/filtered_feature_bc_matrix
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
ge_original_GSC@meta.data <- ge_original_GSC@meta.data %>%
    select(c("orig.ident", "nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"))
rm(ge_original)
gc()

print(paste0("SC1 = ",unique(gte_original$sampleCombined[which(gte_original$neuro.ident == "SC1")])))
print(paste0("SC2 = ",unique(gte_original$sampleCombined[which(gte_original$neuro.ident == "SC2")])))
print(paste0("SC3 = ",unique(gte_original$sampleCombined[which(gte_original$neuro.ident == "SC3")])))
gte_original_GSC <- subset(gte_original, subset = (neuro.ident == "SC1" | neuro.ident == "SC2" | neuro.ident == "SC3"))
gte_original_GSC@meta.data <- gte_original_GSC@meta.data %>%
    select(c("orig.ident", "nCount_RNA", "nFeature_RNA", "log10GenesPerUMI", "mitoRatio"))
rm(gte_original)
gc()

# Create Seurat Objects
gte <- CreateSeuratObject(counts = inter, 
                            project = "Single cell GBM Analysis -- hg38+TE",
                            names.field = 2, names.delim = "-",
                            min.features = 100,
                            min.cells = 1)
gte_combined <- merge(gte, y=gte_original_GSC, 
                        add.cell.ids = c("GBM", "GSC"),
                        project = "Single cell GBM+GSC Analysis -- hg38+TE")
rm(te_mtx)
rm(inter)
rm(gte)
rm(gte_original_GSC)
gc()

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

print("\n")
print("Metadata of combined seurat objects")
gte_combined@meta.data %>% glimpse()
print("\n")
ge_combined@meta.data %>% glimpse()
print("\n")

# Clean metadata
ge_meta <- ge_combined@meta.data
gte_meta <- gte_combined@meta.data

ge_meta <- ge_meta %>% add_column(sampleID = "0")
gte_meta <- gte_meta %>% add_column(sampleID = "0")

# ge_meta$sampleID = ""
# gte_meta$sampleID = ""

print("\n")
print("Separating GBM sample type from cellname")

ge_meta <- ge_meta %>% 
    tibble::rownames_to_column("cells") %>% 
    separate("cells", sep = "[_]", into = c("sample_source", "cells"), remove = FALSE)
gte_meta <- gte_meta %>%  
    tibble::rownames_to_column("cells") %>% 
    separate("cells", sep = "[_]", into = c("sample_source", "cells"), remove = FALSE)

ge_meta %>% glimpse()
gte_meta %>% glimpse()

print("\n")
print("Filling sample IDs")

for (i in 1:length(samples)) {
    length(ge_meta$sampleID[which(ge_meta$sample_source == "GBM" & ge_meta$orig.ident == i)])
    ge_meta$sampleID[which(ge_meta$sample_source == "GBM" & ge_meta$orig.ident == i)] <- samples[i]
    length(gte_meta$sampleID[which(gte_meta$sample_source == "GBM" & gte_meta$orig.ident == i)])
    gte_meta$sampleID[which(gte_meta$sample_source == "GBM" & gte_meta$orig.ident == i)] <- samples[i]
}

print("\n")
print("Filling sampleGSC IDs")

samplesGSC <- samples[(length(samples)-2):length(samples)]
print(samplesGSC)

for (i in 1:length(samplesGSC)) {
    ge_meta$sampleID[which(ge_meta$sample_source == "GSC" & ge_meta$orig.ident == i)] <- samplesGSC[i]
    gte_meta$sampleID[which(gte_meta$sample_source == "GSC" & gte_meta$orig.ident == i)] <- samplesGSC[i]
}

rownames(ge_meta) <- rownames(ge_combined@meta.data)
rownames(gte_meta) <- rownames(gte_combined@meta.data)

ge_combined@meta.data <- ge_meta %>% separate("sampleID", sep = "[_]", into = c("sample_full", "x"), remove = FALSE) %>% select(-x)
gte_combined@meta.data <- gte_meta%>% separate("sampleID", sep = "[_]", into = c("sample_full", "x"), remove = FALSE) %>% select(-x)

## save objects
saveRDS(gte_combined, file= paste0("gte_", args[4]))
print(paste0("Exported initial gte seurat object"))
saveRDS(ge_combined, file=paste0("ge_", args[4]))
print(paste0("Exported initial ge seurat object"))

sessionInfo()
print("Script complete. Exiting")