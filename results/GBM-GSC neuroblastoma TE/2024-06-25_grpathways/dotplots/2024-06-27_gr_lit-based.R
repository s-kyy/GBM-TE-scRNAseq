#!usr/bin/env Rscript

# R version 4.2

## Import 
library(ggplot2)
library(SeuratObject)
library(Seurat)
library(dplyr)
library(gridExtra)
library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr

set.seed(100)
print(getwd())

### Load Data
gbm <- readRDS("D:\\backup files\\2021-11-11\\Cluster\\scratch\\gete-gbm\\results\\GBM-GSC neuroblastoma TE\\2021-09-02\\ge_gbmscIntUmap-subtypes.rds")
DefaultAssay(gbm) <- "integrated"

colnames(gbm@meta.data)
head(gbm$integrated_snn_res.0.3)
head(Idents(gbm))
#setwd("D:\\backup files\\getegbm\\gete-gbm\\results\\GBM-GSC neuroblastoma TE\\2024-06-25_grpathways")

## DEG of resolution 0.3 clusters
Idents(gbm) <- "integrated_snn_res.0.3"

## glutamate associated genes
## Often expressed in glial cells and cancer cells
glu.transporters <- c("SLC1A1", "SLC1A2", "SLC1A3")
grin <- c()
gria <- c()
grm <- c()
grik <- c()

head(genecard.gr$Gene.Symbol)
gr.present <- sort(genecard.gr$Gene.Symbol[ -which(is.na(match(genecard.gr$Gene.Symbol, 
                                                               rownames(gbm)))) ] )
# gr.present <- gr[-which(is.na(match(gr, rownames(gte))))] 
# print(match(gr.present, rownames(gte)))

print(gr.present)

size <- 15

p <- FeaturePlot(gbm, features=gr.present[1:50],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("figs_ge_FeaturePlots_glutamate-receptors_genecard_1-50.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

## save gr.present to csv
write.csv(gr.present, file ="GENECARDS_gr-in-ge.csv", row.names=FALSE)



