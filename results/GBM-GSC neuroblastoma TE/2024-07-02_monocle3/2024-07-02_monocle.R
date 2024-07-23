#!usr/bin/env Rscript

# R version 4.2
set.seed(34)

## Import 
library(ggplot2)
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(100)
print(getwd())

### Load Data
gbm <- readRDS("D:\\backup files\\2021-11-11\\Cluster\\scratch\\gete-gbm\\results\\GBM-GSC neuroblastoma TE\\2021-09-02\\gte_gbmscIntUmap-subtypes.rds")
DefaultAssay(gbm) <- "integrated"
colnames(gbm@meta.data)
head(gbm$integrated_snn_res.0.3)
head(Idents(gbm))
Idents(gbm) <- "integrated_snn_res.0.3"

mono.idents <- read.csv("D:\\backup files\\getegbm\\HMM_CE\\output_gbm\\2023-03-30_sorted_gte_hmm\\predictions.csv")
colnames(mono.idents) <- c("cellnames", "predictions")
gbm@meta.data <- rownames_to_column(gbm@meta.data, "cellnames")
gbm@meta.data <- gbm@meta.data %>% left_join(mono.idents, by = "cellnames")
gbm@meta.data$predictions <- factor(as.character(gbm@meta.data$predictions), levels = sort(unique(gbm@meta.data$predictions)))
Idents(gbm) <- "predictions"

DimPlot(gbm, label=TRUE)

# gbm@meta.data %>% gather(integrated_snn_res.0.3) %>% select(predictions) 
p <- ggplot(gbm@meta.data, aes(x=integrated_snn_res.0.3, y=predictions, fill=predictions)) + 
    geom_bar(position="stack", stat="identity")
    
size <- 5
ggsave("figs_gte_predictions_barplot.tiff", plot = p, units="in", width=size*1.5, height=size*2, dpi=300, compression = 'lzw')

Idents(gbm) <- "integrated_snn_res.0.3"
p <- ggplot(gbm@meta.data, aes(x=predictions, y=integrated_snn_res.0.3, fill=integrated_snn_res.0.3)) + 
    geom_bar(position="stack", stat="identity")
    
size <- 5
ggsave("figs_gte_predictions_barplot_reverse.tiff", plot = p, units="in", width=size*1.5, height=size*2, dpi=300, compression = 'lzw')

Idents(gbm) <- "gbm_subtype"
p <- ggplot(gbm@meta.data, aes(x=predictions, y=gbm_subtype, fill=gbm_subtype)) + 
    geom_bar(position="stack", stat="identity")
    
size <- 5
ggsave("figs_gte_predictions_barplot_gbmsubtype.tiff", plot = p, units="in", width=size*1.5, height=size*2, dpi=300, compression = 'lzw')
