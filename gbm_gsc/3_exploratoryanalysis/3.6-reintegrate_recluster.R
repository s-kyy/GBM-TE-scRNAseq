#!usr/bin/env Rscript

.libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 1 filepath must be supplied: [xxx.rds] [res1] [res2] [sample_name] [fig_dir_name]", call.=FALSE)
} else if (length(args)>=3) {
  
  # verify filepaths
  if (file.exists(args[1])) { 
    obj_path <- args[1] 
    filename <- unlist(strsplit(basename(path_ext_remove(obj_path)), "_"))
    filename <- paste(filename[length(filename)-1],filename[length(filename)], sep="_")
    parent_dir_path_obj <- dirname(obj_path)
    # parent_dir_name_obj <- basename(parent_dir_path_obj)
    res1 <- args[2]
    res2 <- args[3]
  } else {
    stop("Filepaths provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>=4) {
    sample_name <- args[4]
  } else {
    sample_name <- "test"
  }
  if (length(args)>=5) {
    figs_dir_name <- args[5]
  } else {
    figs_dir_name <- "figs_int"
  }

} else {
    stop("Error in calling R or filepaths. Exiting...", call.=FALSE)
}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 

set.seed(108)
options(warn=1) #print warning messages as they occur

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
seurat_obj <- readRDS(obj_path) 
size    <- 5
ndims <- 25

subdir <- file.path(getwd(), paste(format(Sys.Date(), "%Y%m%d"), sample_name, filename, sep="_"))
# subdir <- parent_dir_path_obj

ifelse(!dir.exists(file.path(subdir)),
        dir.create(file.path(subdir),recursive=T),
        "Directory Exists")

figs_dir_path <- file.path(subdir, figs_dir_name)

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

if (grepl("healthy", subdir, fixed = TRUE)) {
  sample_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#CC79A7", "#ff716e",
    "#999999", "#0072B2", "#194c76", 
    "#D55E00", "#3a4f41", "#6699cc", "#713e5a")
  female_samples <- c("SRR9262922", "SRR9262937",
                        "SRR9264382", "SRR9264383",
                        "SRR9264388")

} else if (grepl("gbm", subdir, fixed = TRUE)) {
  sample_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#CC79A7", "#ff716e",
    "#999999", "#0072B2", "#194c76")
  female_samples <- "SF11209"
}

## ========================================= ##
## Preprocess again ##
## ========================================= ##
# Steps: 
# Split by sample, normalize, integrate anew
# Integrated assay : dimensionality reduction and use for UMAP projection and visualization
# RNA assay : Normalize for deg analysis, Scale for heatmap plotting.
# See below for workflow guidance:
# https://github.com/satijalab/seurat/issues/3665
# https://github.com/satijalab/seurat/issues/3775

## Filter genes to those only expressed in 10 or more cells
print(paste("Number of genes in original object:", dim(seurat_obj)[1]))
seurat_obj_new <- CreateSeuratObject(
    counts=seurat_obj$RNA@counts,
    project="GBM_TE-analysis",
    min.cells = 10) # match GBM_Bhaduri
print(paste("Number of genes in filtered object:", dim(seurat_obj)[1]))

#### Compute Quality Metrics & Adjust Metadata ####
meta <- seurat_obj_new@meta.data
meta$cells <- rownames(meta)

# Novelty Score: number of genes per UMI in GE+TE dataset
meta$log10GenesPerUMI <- log10(meta$nFeature_RNA) / log10(meta$nCount_RNA)

# mitochondrial ratio
meta$mitoRatio <- PercentageFeatureSet(object = seurat_obj_new, pattern = "^MT-")
meta$mitoRatio <- meta$mitoRatio / 100

# rename nCount_RNA and nFeature_RNA colnames
meta <- meta %>% 
  dplyr::rename(
    nUMI = nCount_RNA, 
    nGene = nFeature_RNA)

meta$nUMI_prefilt <- seurat_obj$nUMI
meta$nGene_prefilt <- seurat_obj$nGene
meta$log10GenesPerUMI_prefilt <- seurat_obj$log10GenesPerUMI
meta$mitoRatio_prefilt <- seurat_obj$mitoRatio

meta_old <- seurat_obj@meta.data[, -which(names(seurat_obj@meta.data) %in% c("orig.ident", "cells", "nUMI", "nGene", "log10GenesPerUMI", "mitoRatio", "nCount_RNA", "nFeature_RNA", "seurat_clusters"))]

# Transfer metadata from old object
seurat_obj_new@meta.data <- cbind(meta, meta_old)
rm("meta_old", "seurat_obj", "meta")
seurat_obj <- seurat_obj_new
rm("seurat_obj_new")

## Split by Sample
DefaultAssay(seurat_obj) <- "RNA"
print(paste("unique samples to split by:", unique(seurat_obj$sample)))
obj_list <- SplitObject(seurat_obj, split.by="sample")
print(paste("Length of list", length(obj_list)))
rm("seurat_obj")

#### ===================================================================== ####
#### Normalize the RNA reads per tumor sample via Log Normalization method ####
#### ===================================================================== ####
p_list <- vector(mode="list",length(obj_list))
names(p_list) <- sort(names(obj_list))

for(i in 1:length(obj_list)){
    print(paste0("Processing Sample: ", obj_list[[i]]$sample[1]))
    DefaultAssay(obj_list[[i]]) <- "RNA"
    # default: LogNormalize, save to : $RNA@data
    obj_list[[i]] <- NormalizeData(obj_list[[i]], verbose = TRUE) 
    obj_list[[i]] <- FindVariableFeatures(obj_list[[i]], 
                              selection.method="vst",
                              nfeatures=2500, #default=2000
                              verbose = TRUE)
    #Create Variable Feature Plots
    top10 <- head(VariableFeatures(obj_list[[i]]), 10)
    p <- VariableFeaturePlot(obj_list[[i]])
    p <- LabelPoints(plot = p, points = top10, repel = TRUE)
    ggsave(file.path(figs_dir_path, paste0(filename,obj_list[[i]]$sample[1],"_varfeatgenes.tiff")), 
          plot = p, units="in", width=size*1.1, height=size*0.8, dpi=300, compression = 'lzw')
    p_list[[obj_list[[i]]$sample[1]]] <- p
}

saveRDS(p_list, file=file.path(figs_dir_path,paste0(filename,"_allvarfeatplots.rds")))
print("Variable Features computed")

#### =========================================== ####
#### Integrate Seurat Objects together: 30min-1h30min ####
#### =========================================== ####
## Anchors across all samples will be used to integrated each sample
## This step helps reduce variation introduced across donor samples 
## (i.e. variation in tumour extraction, experimental processing of samples, 
## gender-based differences and so on)
## Stuck with CCA method to find anchors, because the datasets are a manageable size (<50K cells)
obj_anchors <- FindIntegrationAnchors(
                  object.list = obj_list, 
                  dims=1:ndims,
                  normalization.method = "LogNormalize", 
                  anchor.features = 2500, #match number in FindVariableFeatures
                  )
seurat_obj <- IntegrateData(anchorset = obj_anchors, dims=1:ndims)
# saveRDS(seurat_obj, file = file.path(subdir, paste0(filename, "_int.rds")) )
print("Integration Complete")

## ========================================= ##
## Scale and Reduce Dimensionality without regression ##
## ========================================= ##
DefaultAssay(seurat_obj) <- "integrated"
# Scale Data ~5min, <200Mb
seurat_obj <- ScaleData(object = seurat_obj, verbose = FALSE) 

seurat_obj@meta.data$sample <- factor(seurat_obj@meta.data$sample, 
                                        levels = sort(unique(seurat_obj@meta.data$sample)))
# Run PCA and UMAP ~20min
seurat_obj <- RunPCA(seurat_obj, 
                      npcs = ndims, 
                      verbose = FALSE)

#### =========================================== ####
#### Test PCA levels ####
#### =========================================== ####
# Determine percent of variation associated with each PC
DefaultAssay(seurat_obj) <- "integrated"
pct_var_per_pc <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cum_pct_per_pc <- cumsum(pct_var_per_pc)

# Determine which PC exhibits a cumulative percentage of variation 
# greater than 90% and variation associated with the PC is less than 5%
pc_most_var <- which(cum_pct_per_pc > 90 & pct_var_per_pc < 5)[1]
print(paste("Minimum PC that retains more than 90% variation and less than 5% variation compared to the next PC:", pc_most_var))

# Determine the difference between variation of PC and subsequent PC
# last point where change of % of variation is more than 0.1%.
pc_10 <- sort(which((pct_var_per_pc[1:length(pct_var_per_pc) - 1] - pct_var_per_pc[2:length(pct_var_per_pc)]) > 0.1), decreasing = T)[1] + 1
print(paste("Minimum PC with a difference in variation of 0.1% compared to next PC:", pc_10))

# Minimum of the two calculation
min_pc <- min(pc_most_var, pc_10)
print(paste("Minumum PC between the two options:", min_pc))

plot_df <- data.frame(dimensions = 1:length(pct_var_per_pc),
           stdev = seurat_obj[["pca"]]@stdev,
           pct_var_per_pc = pct_var_per_pc,
           cum_pct_per_pc = cum_pct_per_pc)
print(plot_df)
write.csv(plot_df, file.path(figs_dir_path, paste0(filename, "_pca_filtDf.csv") ))

# Plot % variation to Elbowplot (modified from Seurat Elbow Plot Function)
p <- plot_df %>% ggplot(aes(x = dimensions, y = stdev)) +
    geom_point() +
    labs(x = "", y = "Standard Deviation") +
    geom_text(
      label=format(round(cum_pct_per_pc, 1), nsmall = 1), 
      nudge_x = 0.5, nudge_y = 0.5, 
      check_overlap = T,
      size=2) +
    theme_classic() 
ggsave(file.path(figs_dir_path, paste0(filename,"_elbow_filtDf.tiff")), 
  plot = p, units="in", width=size*0.7, height=size*0.7, dpi=300, compression = 'lzw')
print("Exported ElbowPlot")

#### =========================================== ####
#### Output UMAP plots by sample ####
#### =========================================== ####

### Project Minimum number of PCs to UMAP plot
seurat_obj <- RunUMAP(seurat_obj, 
                          reduction = "pca", 
                          dims = 1:min_pc, 
                          umap.method = "uwot", 
                          metric = "cosine")
# saveRDS(seurat_obj, file = file.path(subdir,paste0(filename, "_filtDf_umap.rds")) )
# print(paste("ScaleData, RunPCA, and RunUMAP completed in:", subdir))
print(paste("RunUMAP completed"))

p <- DimPlot(seurat_obj, reduction = "umap", group.by = "sample") 
ggsave(file.path(figs_dir_path, paste0(filename,"_UMAP_sample_filtDoublets.tiff")), 
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')
print("Exported UMAP by sample")


#### =========================================== ####
#### Cluster: K-nearest neighbor graph
#### =========================================== ####
# # rename old cluster annotations
# metadata <- seurat_obj@meta.data
# metadata <- rename(metadata, filtDF_int_snn_res03 = "integrated_snn_res.0.3")
# metadata <- rename(metadata, filtDF_int_snn_res04 = "integrated_snn_res.0.4")
# seurat_obj@meta.data <- metadata
# glimpse(metadata)
# rm("metadata")

## identify clusters
seurat_obj <- FindNeighbors(seurat_obj,dims=1:min_pc,reduction="pca")
seurat_obj <- FindClusters(seurat_obj, resolution = res1)
seurat_obj <- FindClusters(seurat_obj, resolution = res2)
cat("KNN clustering complete. Saved seurat objects to", subdir,"\n")

#### =========================================== ####
#### Find unique markers per cluster ####
#### =========================================== ####
## Use normalized counts for DGE (Luecken & Theis 2019)
DefaultAssay(seurat_obj) <- "RNA" 

# Use normalized counts to perform differential gene analysis
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE) # default: LogNormalize
saveRDS(seurat_obj, file = file.path(subdir, paste0(filename, "_int.rds")))

seurat_obj <- ScaleData(seurat_obj,verbose = FALSE)
saveRDS(GetAssayData(seurat_obj,slot='scale.data'), file = file.path(subdir, paste0(filename, "_scaledata.rds")))
print("Normalization and Centering + Scaling Complete")

print("Finding Markers...")
Idents(seurat_obj) <- paste0("integrated_snn_res.", res1)
markers <- FindAllMarkers(
  seurat_obj, 
  only.pos = FALSE, 
  min.pct = 0.25, 
  logfc.threshold = 0.1)
write.csv(markers, file = file.path(subdir, paste0(filename, "_markers_",res1,".csv")), row.names=TRUE)
rm("markers")

Idents(seurat_obj) <- paste0("integrated_snn_res.", res2)
markers <- FindAllMarkers(
  seurat_obj, 
  only.pos = FALSE, 
  min.pct = 0.25, 
  logfc.threshold = 0.1)
write.csv(markers, file = file.path(subdir, paste0(filename, "_markers_",res2,".csv")), row.names=TRUE)
rm("markers")

#### End of Script
sessionInfo()