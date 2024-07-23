#!usr/bin/env Rscript

#### Import Packages ####
set.seed(34)

library(Seurat)
library(Matrix)
library(ggplot2)
library(genefilter)
library(tidyverse)
library(scales)
library(AnnotationHub)

set.seed(34)

#### Load Datasets ####
gte <- readRDS("temp/gte.rds")
ge <- readRDS("temp/ge.rds")

#### Quality Control ####
filt_ge <- subset(x = ge, subset= (nUMI >= 1000) & 
                            (nGene >= 300) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

barcode <- colnames(filt_ge)
filt_gte <- subset(gte, cells=barcode)

cat(paste0("\n ge contains", 
           dim(filt_ge), " genes x cells\n"))
cat(paste0("Filtered out ", 
           dim(ge)[2] - dim(filt_ge)[2], " cells in ge\n" ))

cat(paste0("\n gte contains", 
           dim(filt_gte), " genes x cells\n"))
cat(paste0("Filtered out ", 
           dim(gte)[2] - dim(filt_gte)[2], " cells in ge\n" ))

saveRDS(filt_gte, file = "temp/gte.qc.rds")
saveRDS(filt_ge,  file = "temp/ge.qc.rds")

rm("ge")
rm("gte")
gc()

#### Normalize each dataset individually ####


#### Integrate Seurat Objects together ####


#### End of Script ####
sessionInfo()
