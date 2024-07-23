#!usr/bin/env Rscript

# INPUTS
# - samples.csv = csv file with 1 column of strings 
# containing desired samples and no headers  
# - feature_bc_matrix = One mapped to reference genes, 
# Second mapped to transposable element (TE)/retrotransposon (RT) annotations

#### Import Packages ####
set.seed(34)

library(Seurat)
library(Matrix)
library(genefilter)
library(tidyverse)
library(scales)
library(AnnotationHub)

set.seed(34)

#### Load Datasets ####
sampleNames <- read.csv("./samples.csv", header = FALSE)  
sampleNames <- sampleNames[,1]
matrix.path <- "./aggr/outs/count/feature_bc_matrix"
matrix.path.RT <- "./aggr.rt/outs/count/feature_bc_matrix"

reads <- Read10X(matrix.path)
reads.RT <- Read10X(matrix.path.RT)
cat(dim(reads))
cat(dim(reads.RT))

#### Create Seurat Objects ####

# find intersecting cell names between genes and retrotransposon datasets 
cellnames <- intersect(colnames(reads),colnames(reads.RT))
length(cellnames) 

# Merge sparse matrices with intersecting cell names
intersect <- rbind(reads[,cellnames],reads.RT[,cellnames])
rm("reads.RT")
gc()

gte <- CreateSeuratObject(
    counts=intersect,
    project="reference-genes_retrotransposons", 
    names.field = 2,
    names.delim = "-", 
    min.features = 100, 
    min.cells = 1)
# gte.orig_id <- gte@meta.data$orig.ident
# names(gte.orig_id) <- rownames(gte@meta.data)
rm("intersect")
gc()

ge <- CreateSeuratObject(
    counts=reads,
    project="reference-genes",
    names.field = 2,
    names.delim = "-", 
    min.features = 100, 
    min.cells = 1)
# ge.orig_id <- ge@meta.data$orig.ident
# names(ge.orig_id) <- rownames(ge@meta.data)
rm("reads")
gc()

#### Create Results Folder with today's date ####
cat("Making Today's directory\n")
cat(paste0("Current working directory: ", getwd(),"\n"))

subdir <- format(Sys.Date(), "%Y%m%d")
dir.create(file.path(getwd(), subdir))
setwd(file.path(getwd(), subdir))
print(paste0("New working directory: ", getwd()))

dir.create(file.path(getwd(), "temp"))
saveRDS(gte, file = "temp/gte.rds")
saveRDS(ge,  file = "temp/ge.rds")

barcode <- colnames(gte)
ge <- subset(ge, cells=barcode)
cat(paste0("gte genes x cellbarcodes: ", dim(gte), "\n"))
cat(paste0("ge genes x cellbarcodes: ", dim(ge), "\n"))

rm("barcode")
gc()

#### Compute Quality Metrics & Adjust Metadata ####

# make tmp metadata variables
metadata_gte <- gte@meta.data
metadata_ge <- ge@meta.data

# rename nCount_RNA and nFeature_RNA colnames
metadata_gte <- metadata_gte %>% dplyr::rename(nUMI = nCount_RNA,
                                               nGene = nFeature_RNA)
metadata_ge <- metadata_ge %>% dplyr::rename(nUMI = nCount_RNA,
                                             nGene = nFeature_RNA)

### Add cell IDs (rownames) to metadata
metadata_gte$cells <- rownames(metadata_gte)
metadata_ge$cells <- rownames(metadata_ge)

# number of genes per UMI in GE+TE dataset
metadata_gte$log10GenesPerUMI <- log10(metadata_gte$nFeature_RNA) / log10(metadata_gte$nCount_RNA)
metadata_ge$log10GenesPerUMI <- log10(metadata_ge$nFeature_RNA) / log10(metadata_ge$nCount_RNA)

# mitochondrial ratio
metadata_gte$mitoRatio <- PercentageFeatureSet(object = gte, pattern = "^MT-")
metadata_gte$mitoRatio <- metadata_gte$mitoRatio / 100
metadata_ge$mitoRatio <- PercentageFeatureSet(object = ge, pattern = "^MT-")
metadata_ge$mitoRatio <- metadata_ge@meta.data$mitoRatio / 100

### Fill "sample" column
metadata_gte$sample <- NA
metadata_ge$sample <- NA

for(i in 1:length(sampleNames)) {
    metadata_gte$sample[which(metadata_gte$orig.ident == i)] <- sampleNames[i]
    metadata_ge$sample[which(metadata_ge$orig.ident == i)] <- sampleNames[i]
}


head(metadata_ge)
head(metadata_gte)

gte@meta.data <- metadata_gte
ge@meta.data <- metadata_ge

saveRDS(gte, file = "temp/gte.rds")
saveRDS(ge,  file = "temp/ge.rds")

#### End of Script ####
sessionInfo()
