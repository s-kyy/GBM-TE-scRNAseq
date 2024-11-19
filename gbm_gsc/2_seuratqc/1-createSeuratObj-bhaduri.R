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
library(tidyverse)

set.seed(34)

#### Load Datasets ####
sample_names <- read.csv(file.path("..","0_downloads","2021-03-02_bhaduri","samples.csv"), header = false) # USER ADJUST
sample_names <- sample_names[,1]
matrix_path_GE <- file.path("..","1_scrna-seq_mapping","2021-03-10_bhaduri_aggr_ge","aggr","outs","count","feature_bc_matrix") # USER ADJUST
matrix_path_TE <- file.path("..","1_scrna-seq_mapping","2021-03-25_bhaduri_aggr_te","aggr","outs","count","feature_bc_matrix") # USER ADJUST

reads_GE <- Read10X(matrix_path_GE)
reads_TE <- Read10X(matrix_path_TE)
cat(dim(reads_GE))
cat(dim(reads_TE))
rm("matrix_path_GE")
rm("matrix_path_TE")

print(paste0("Total number of cells (reads_GE) before qc filtering: ", 
              length(reads_GE[1,]),"\n"))
              # 60065 cells in the original TE-mapped sparse matrix
print(paste0("Total number of cells (reads_TE) before qc filtering: ", 
              length(reads_TE[1,]),"\n"))
              # 37204 cells in the original TE-mapped sparse matrix

# loop through each colname (each cellbarcode) in reads_TE matrix to match those of normal reads matrix
for(i in 1:length(reads_TE[1,])){
    if ( str_detect(colnames(reads_TE)[i], "-1$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-1$", "-11")
    } else if ( str_detect(colnames(reads_TE)[i], "-2$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-2$", "-5")
    } else if ( str_detect(colnames(reads_TE)[i], "-3$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-3$", "-1")
    } else if ( str_detect(colnames(reads_TE)[i], "-4$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-4$", "-12")
    } else if ( str_detect(colnames(reads_TE)[i], "-5$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-5$", "-6")
    } else if ( str_detect(colnames(reads_TE)[i], "-6$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-6$", "-2")
    } else if ( str_detect(colnames(reads_TE)[i], "-7$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-7$", "-3")
    } else if ( str_detect(colnames(reads_TE)[i], "-8$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-8$", "-9")
    } else if ( str_detect(colnames(reads_TE)[i], "-9$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-9$", "-8")
    } else if ( str_detect(colnames(reads_TE)[i], "-10$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-10$", "-7")
    } else if ( str_detect(colnames(reads_TE)[i], "-11$") ) {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-11$", "-4")
    } else {
        colnames(reads_TE)[i] <- str_replace(colnames(reads_TE)[i], "-12$", "-10")        
    }
}

#### Create Seurat Objects ####

# find intersecting cell names between genes and retrotransposon datasets 
cellnames <- intersect(colnames(reads_GE),colnames(reads_TE))
print(paste0("Total number of cells (reads_TE) before qc filtering: ", 
              length(cellnames)))
              # 37104 cells will be present in the combined sparse matrix

# Merge sparse matrices with intersecting cell names
intersect <- rbind(reads_GE[,cellnames],reads_TE[,cellnames])
rm("reads_TE")
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
    counts=reads_GE,
    project="reference-genes",
    names.field = 2,
    names.delim = "-", 
    min.features = 100, 
    min.cells = 1)
# ge.orig_id <- ge@meta.data$orig.ident
# names(ge.orig_id) <- rownames(ge@meta.data)
rm("reads_GE")
gc()

#### Create Results Folder with today's date ####
cat("Making Today's directory\n")
cat(paste0("Current working directory: ", getwd(),"\n"))

subdir <- paste0(format(Sys.Date(), "%Y%m%d"), "_", "bhaduriGBM")
ifelse(!dir.exists(file.path(getwd(),subdir, "temp")),
        dir.create(file.path(getwd(),subdir, "temp"),recursive=T),
        "Directory Exists")
saveRDS(gte, file = file.path(getwd(),subdir,"temp","gte.rds"))
saveRDS(ge, file = file.path(getwd(),subdir,"temp","ge.rds"))

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

### Fill "sample_orig" column - based on NCBI page
metadata_gte$sample_orig <- NA
metadata_ge$sample_orig <- NA
metadata_gte$sample <- NA
metadata_ge$sample <- NA
 
index <- (2, 3, 4, 12, 6, 8, 11, 5, 1, 10, 7, 9)
sample_names <- sample_names[index]

for(i in 1:length(sample_names)) {
    metadata_gte$sample__orig[which(metadata_gte$orig.ident == i)] <- sample_names[i]
    metadata_ge$sample__orig[which(metadata_ge$orig.ident == i)] <- sample_names[i]
}

### Fill "sampleCombined" column - based on patient 
for(i in 1:length(metadata_gte$sample_orig)) {
    if ( str_detect(metadata_gte$sample__orig[i], "_1$") ) {
        metadata_ge$sample[i] <- str_replace(metadata_gte$sample__orig[i], "_1$", "")
    } else {
        metadata_ge$sample[i] <- str_replace(metadata_gte$sample__orig[i], "_2$", "")
}

head(metadata_ge)
head(metadata_gte)

gte@meta.data <- metadata_gte
ge@meta.data <- metadata_ge

saveRDS(gte, file = file.path(getwd(),subdir,"temp","gte.rds"))
saveRDS(ge,  file = file.path(getwd(),subdir,"temp","ge.rds"))
ifelse(!dir.exists(file.path(getwd(),subdir, "seurat_obj")),
        dir.create(file.path(getwd(),subdir, "seurat_obj"),recursive=T),
        "Directory Exists")
file.rename(from = file.path(getwd(),subdir,"temp","gte.rds"), 
            to = file.path(getwd(),subdir,"seurat_obj","gte.rds"))
file.rename(from = file.path(getwd(),subdir,"temp","ge.rds"), 
            to = file.path(getwd(),subdir,"seurat_obj","ge.rds"))

#### End of Script ####
sessionInfo()
