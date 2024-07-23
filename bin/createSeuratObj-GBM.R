#!usr/bin/env Rscript

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
gbm.matrix.path <- ""
gbm.matrix.path.RT <- ""
gsc.matrix.path <- ""
gsc.matrix.path.RT <- ""
# e.g. "~/.../aggr/outs/count/filtered_feature_bc_matrix"

gbm.reads <- Read10X(gbm.matrix.path)
gbm.reads.RT <- Read10X(gbm.matrix.path.RT)
cat(dim(gbm.reads))
cat(dim(gbm.reads.RT))

gsc.reads <- Read10X(gsc.matrix.path)
gsc.reads.RT <- Read10X(gsc.matrix.path.RT)
cat(dim(gsc.reads))
cat(dim(gsc.reads.RT))

#### Create Seurat Objects ####

# GBM: find intersecting cell names between genes and retrotransposon datasets 
cellnames <- intersect(colnames(gbm.reads),colnames(gbm.reads.RT))
length(cellnames) 

# Merge sparse matrices with intersecting cell names
gbm.intersect <- rbind(gbm.reads[,cellnames],gbm.reads.RT[,cellnames])
rm(gbm.reads.RT)
gc()

gbm.gte <- CreateSeuratObject(
    counts=gbm.intersect,
    project="gbm_GenesRetrotransposons",
    names.field = 2,
    names.delim = "-", 
    min.features = 100, 
    min.cells = 1)
# gbm.gte.orig_id <- gbm.gte@meta.data$orig.ident
# names(gbm.gte.orig_id) <- rownames(gbm.gte@meta.data)
rm(gbm.intersect)
gc()

gbm.ge <- CreateSeuratObject(
    counts=gbm.reads,
    project="gbm_Refgenes",
    names.field = 2,names.delim = "-", min.features = 100, min.cells = 1)
# gbm.ge.orig_id <- gbm.ge@meta.data$orig.ident
# names(gbm.ge.orig_id) <- rownames(gbm.ge@meta.data)
rm(gbm.reads)
gc()

# GSC: find intersecting cell names between genes and retrotransposon datasets
cellnames <- intersect(colnames(gsc.reads),colnames(gsc.reads.RT))
length(cellnames) 

# Merge sparse matrices with intersecting cell names
gsc.intersect <- rbind(gsc.reads[,cellnames],gsc.reads.RT[,cellnames])
rm(gsc.reads.RT)
gc()

gsc.gte <- CreateSeuratObject(
    counts=gsc.intersect,
    project="neuroblastoma_GenesRetrotransposons",
    names.field = 2,
    names.delim = "-", 
    min.features = 100, 
    min.cells = 1)
# gsc.gte.orig_id <- gsc.gte@meta.data$orig.ident
# names(gsc.gte.orig_id) <- rownames(gsc.gte@meta.data)
rm(gsc.intersect)
gc()

gsc.ge <- CreateSeuratObject(
    counts=gsc.reads,
    project="neuroblastoma_Refgenes",
    names.field = 2,names.delim = "-", min.features = 100, min.cells = 1)
# gbm.ge.orig_id <- gbm.ge@meta.data$orig.ident
# names(gbm.ge.orig_id) <- rownames(gbm.ge@meta.data)
rm(gsc.reads)
gc()
    
#### Create Results Folder with today's date ####
cat("Making Today's directory\n")
cat(paste0("Current working directory: ", getwd(),"\n"))

subdir <- format(Sys.Date(), "%Y%m%d")
dir.create(file.path(getwd(), subdir))
setwd(file.path(getwd(), subdir))
print(paste0("New working directory: ", getwd()))

dir.create(file.path(getwd(), "temp"))
saveRDS(gbm.gte, file = "temp/gbm.gte.rds")
saveRDS(gbm.ge,  file = "temp/gbm.ge.rds")
saveRDS(gsc.gte, file = "temp/gsc.gte.rds")
saveRDS(gsc.ge,  file = "temp/gsc.ge.rds")

barcode <- colnames(gbm.gte)
gbm.ge <- subset(gbm.ge, cells=barcode)
cat(paste0("GBM.gte genes x cellbarcodes: ", dim(gbm.gte), "\n"))
cat(paste0("GBM.ge genes x cellbarcodes: ", dim(gbm.ge), "\n"))

barcode <- colnames(gsc.gte)
gsc.ge <- subset(gsc.ge, cells=barcode)
cat(paste0("GSC.gte genes x cellbarcodes: ", dim(gsc.gte), "\n"))
cat(paste0("GSC.ge genes x cellbarcodes: ", dim(gsc.ge), "\n"))



#### Normalize each dataset individually ####


#### Integrate Seurat Objects together ####


