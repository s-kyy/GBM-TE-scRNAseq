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
filt.ge <- subset(x = ge, subset= (nUMI >= 1000) & 
                            (nGene >= 300) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

barcode <- colnames(filt.ge)
filt.gte <- subset(gte, cells=barcode)

cat(paste0("\n ge contains", 
           dim(filt.ge), " genes x cells\n"))
cat(paste0("Filtered out ", 
           dim(ge)[2] - dim(filt.ge)[2], " cells in ge\n" ))

cat(paste0("\n gte contains", 
           dim(filt.gte), " genes x cells\n"))
cat(paste0("Filtered out ", 
           dim(gte)[2] - dim(filt.gte)[2], " cells in ge\n" ))

saveRDS(filt.gte, file = "temp/gte.qc.rds")
saveRDS(filt.ge,  file = "temp/ge.qc.rds")
cat("Filtered & saved seurat objects to temp/\n")

rm("ge")
rm("gte")
gc()

#### Create Unique Identities for each Sample ####
# needed if merging from multiple cellranger aggr outputs
# todo

#### Normalize each dataset individually ####

gte.temp <- CreateAssayObject(counts = filt.gte@assays[["RNA"]]@counts)

gte.temp <- NormalizeData(gte.temp, 
                          verbose = T, 
                          normalization.method = "LogNormalize") %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData()

filt.gte@assays[["RNA"]]@counts <- gte.temp@counts
filt.gte@assays[["RNA"]]@scale.data <- gte.temp@scale.data

ge.temp  <- CreateAssayObject(counts = filt.ge@assays[["RNA"]]@counts)
ge.temp <- NormalizeData(ge.temp, 
                          verbose = T, 
                          normalization.method = "LogNormalize") %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData()

filt.ge@assays[["RNA"]]@counts <- ge.temp@counts
filt.ge@assays[["RNA"]]@scale.data <- ge.temp@scale.data

rm("ge.temp")
rm("gte.temp")
gc()

gte.list <- SplitObject(object=filt.gte, split.by="orig.ident")
ge.list <- SplitObject(object=filt.ge, split.by="orig.ident")

rm("filt.gte")
rm("filt.ge")
gc()


# Change active.ident as orig.ident for GBM samples
for (i in 1:length(gte.list)) {
    gte.list[[i]]@active.ident <- factor(gte.list[[i]]@meta.data$orig.ident)
    ge.list[[i]]@active.ident <- factor(ge.list[[i]]@meta.data$orig.ident)
}

# Normalize the RNA reads per tumor sample via 
# Log Normalization method, then remove most variable features (outliers)

for(i in 1:length(ge.list)){
    
    print(paste0("Processing Sample Number: ",i))
    
    DefaultAssay(gte.list[[i]]) <- "RNA"
    DefaultAssay(ge.list[[i]]) <- "RNA"
    
    gte.list[[i]] <- NormalizeData(gte.list[[i]],verbose=TRUE)
    ge.list[[i]]  <- NormalizeData(ge.list[[i]],verbose=TRUE)
    
    gte.list[[i]] <- FindVariableFeatures(gte.list[[i]],
                                          selection.method="vst",
                                          nfeatures=2500,
                                          verbose=TRUE)
    ge.list[[i]]  <- FindVariableFeatures(ge.list[[i]],
                                          selection.method="vst",
                                          nfeatures=2500,
                                          verbose=TRUE)
}
cat("Most Variable Features computed\n")

#### Integrate Seurat Objects together (memory: 60Gb) ~1h30min ####
gte.anchors <- FindIntegrationAnchors(gte.list,dims=1:20)
gte.integrated <- IntegrateData(anchorset = gte.anchor,dims=1:20)

ge.anchors <- FindIntegrationAnchors(ge.list,dims=1:20) 
ge.integrated <- IntegrateData(anchorset = ge.anchor,dims=1:20)
cat("Most Integration Complete")

saveRDS(gte.integrated, "temp/gte.integrated.rds")
saveRDS(ge.integrated, "temp/ge.integrated.rds")
cat("Saved integrated seurat objects to temp/")

rm("gte.list")
rm("ge.list")
rm("gte.anchors")
rm("ge.anchors")
gc()

#### Create Variable Feature Plots ####
DefaultAssay(object = gte.integrated) <- "integrated"
DefaultAssay(object = ge.integrated) <- "integrated"

# Identify the 10 most highly variable genes
size    <- 7
top10   <- head(VariableFeatures(gte.integrated), 10)
p       <- VariableFeaturePlot(gte.integrated)
p2      <- LabelPoints(plot = p, points = top10, repel = TRUE)
ggsave("temp/gte_variablegenes.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

top10   <- head(VariableFeatures(ge.integrated), 10)
p       <- VariableFeaturePlot(ge.integrated)
p2      <- LabelPoints(plot = p, points = top10, repel = TRUE)
ggsave("temp/ge_variablegenes.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

#### Scale, PCA, UMAP ####
# Standard workflow for visualization and clustering 
# (ScaleData, RunPCA, RunUMAP, FindNeighbors, FindClusters)
gte.integrated <- ScaleData(object = gte.integrated, verbose = FALSE)
ge.integrated <- ScaleData(object = ge.integrated, verbose = FALSE)

gte.integrated@meta.data$sample <- factor(gte.integrated@meta.data$sample, 
                                          levels = unique(gte.integrated@meta.data$sample))
ge.integrated@meta.data$sample <- factor(ge.integrated@meta.data$sample, 
                                          levels = unique(ge.integrated@meta.data$sample))

# Run PCA and UMAP ~20min
gte.integrated.pca <- RunPCA(object = gte.integrated, npcs = 20, verbose = FALSE)
ge.integrated.pca  <- RunPCA(object = ge.integrated, npcs = 20, verbose = FALSE)

rm("gte.integrated")
rm("ge.integrated")
gc()

gte.integrated.pca <- RunUMAP(object = gte.integrated.pca, 
                              reduction = "pca", 
                              dims = 1:20, 
                              umap.method = "uwot", 
                              metric = "cosine")

ge.integrated.pca <- RunUMAP(object = ge.integrated.pca, 
                              reduction = "pca", 
                              dims = 1:20, 
                              umap.method = "uwot", 
                              metric = "cosine")

saveRDS(gte.integrated.pca, "temp/gte.integrated.umap.rds")
saveRDS(ge.integrated.pca, "temp/ge.integrated.umap.rds")
cat("ScaleData, RunPCA, and RunUMAP complete.\n Saved Seurat objects to temp/")

# Determine the K-nearest neighbor graph (with first 20 PCAs)
gte.integrated.pca <- FindNeighbors(gte.integrated.pca,dims=1:20,reduction="pca")
gte.integrated.pca <- FindClusters(gte.integrated.pca, resolution = 0.3)
gte.integrated.pca <- FindClusters(gte.integrated.pca, resolution = 0.4)
saveRDS(gte.integrated.pca, "temp/gte.integrated.umap.clustered.rds")
cat("KNN clustering complete. Saved seurat objects to temp/")

ge.integrated.pca <- FindNeighbors(ge.integrated.pca,dims=1:20,reduction="pca")
ge.integrated.pca <- FindClusters(ge.integrated.pca, resolution = 0.3)
ge.integrated.pca <- FindClusters(ge.integrated.pca, resolution = 0.4)
saveRDS(ge.integrated.pca, "temp/ge.integrated.umap.clustered.rds")
cat("KNN clustering complete. Saved seurat objects to temp/")

size <- 5
p <- DimPlot(gte.integrated.pca, reduction = "umap", group.by = "sample") + 
    ggtitle("GRCh38+TE by Sample Origin")
ggsave("temp/gte_UMAP-sample.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

p <- DimPlot(ge.integrated.pca, reduction = "umap", group.by = "sample") + 
    ggtitle("GRCh38 by Sample Origin")
ggsave("temp/ge_UMAP-sample.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

p <- DimPlot(gte.integrated.pca, reduction = "umap", group.by = "cluster") + 
    ggtitle("GRCh38+TE by Sample Origin")
ggsave("temp/gte_UMAP-cluster.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

p <- DimPlot(ge.integrated.pca, reduction = "umap", group.by = "cluster") + 
    ggtitle("GRCh38 by Sample Origin")
ggsave("temp/ge_UMAP-cluster.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

cat("Seurat object processing complete. Saved to temp/")
#### End of Script ####
sessionInfo()
