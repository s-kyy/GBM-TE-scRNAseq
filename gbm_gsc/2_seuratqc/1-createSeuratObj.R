#!usr/bin/env Rscript

# INPUTS
# - samples.csv = csv file with 1 column of strings 
# containing desired samples and no headers  
# - feature_bc_matrix = One mapped to reference genes, 
# Second mapped to transposable element (TE)/retrotransposon (RT) annotations

#### Parse Arguments ####
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<=2) {
  stop("At least 3 filepaths must be supplied for [sample.csv], [GE-feature_bc_matrix], [TE-feature_bc_matrix]. Optional args include: [output-dir-name] [output-prefix]", call.=FALSE)
} else if (length(args)>=3) {
  # verify paths
  if (file.exists(args[1])){ sample_path <- args[1]  } 
  if (dir.exists(args[2])){ matrix_path_GE_parent <- args[2]  } 
  if (dir.exists(args[3])){ matrix_path_TE_parent <- args[3]  }
  output_dir_name <- ""
  output_file_GTE <- "gte.rds"  
  output_file_GE <- "ge.rds"  
  if (length(args)==4) { 
    output_dir_name <- args[4]
  } else if (length(args)==5) { 
    output_dir_name <- args[4]
    output_file_GTE <- paste0(args[5], "_", "gte.rds") 
    output_file_GE <- paste0(args[5], "_", "ge.rds") 
  }
}

#### Import Packages ####
set.seed(34)

library(Seurat)
library(Matrix)
library(tidyverse)

set.seed(34)

#### Load Datasets ####
sampleNames <- read.csv(sample_path, header = FALSE) 
sampleNames <- sampleNames[,1]
matrix_path_GE <- file.path( matrix_path_GE_parent, "aggr","outs","count","feature_bc_matrix") 
matrix_path_TE <- file.path( matrix_path_TE_parent, "aggr","outs","count","feature_bc_matrix") 

reads_GE <- Read10X(matrix_path_GE)
reads_TE <- Read10X(matrix_path_TE)
cat(dim(reads_GE))
cat(dim(reads_TE))

print(paste0("Total number of cells (reads_GE) before qc filtering: ", 
              length(reads_GE[1,]),"\n"))
print(paste0("Total number of cells (reads_TE) before qc filtering: ", 
              length(reads_TE[1,]),"\n"))

#### Create Seurat Objects ####

# find intersecting cell names between genes and retrotransposon datasets 
cellnames <- intersect(colnames(reads_GE),colnames(reads_TE))
length(cellnames) 

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

subdir <- paste0(format(Sys.Date(), "%Y%m%d"), "_", output_dir_name)
ifelse(!dir.exists(file.path(getwd(),subdir, "temp")),
        dir.create(file.path(getwd(),subdir, "temp"),recursive=T),
        "Directory Exists")
saveRDS(gte, file = file.path(getwd(),subdir,"temp",output_file_GTE))
saveRDS(ge, file = file.path(getwd(),subdir,"temp",output_file_GE))

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

saveRDS(gte, file = file.path(getwd(),subdir,"temp",output_file_GTE))
saveRDS(ge,  file = file.path(getwd(),subdir,"temp",output_file_GE))
file.rename(from = file.path(getwd(),subdir,"temp",output_file_GTE), 
            to = file.path(getwd(),subdir,output_file_GTE))
file.rename(from = file.path(getwd(),subdir,"temp",output_file_GE), 
            to = file.path(getwd(),subdir,output_file_GE))

#### End of Script ####
sessionInfo()
