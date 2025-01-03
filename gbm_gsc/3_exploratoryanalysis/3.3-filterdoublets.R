#!usr/bin/env Rscript

.libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 1 filepath must be supplied: [xxx.rds] [yyy.rds] [sample_name]", call.=FALSE)
} else if (length(args)>=2) {
  
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2])) { 
    obj_path <- args[1] 
    obj_path_2 <- args[2] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(dirname(obj_path))
    parent_dir_name_obj <- basename(parent_dir_path_obj)
  } else {
    stop("Filepaths provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>1) {
    sample_name <- args[3]
  } else {
    sample_name <- ""
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

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
seurat_obj <- readRDS(obj_path) 
seurat_obj_2 <- readRDS(obj_path_2)
size    <- 5
ndims <- 25

subdir <- file.path(getwd(), paste0(format(Sys.Date(), "%Y%m%d"), "_", sample_name,"_", filename))
# subdir <- parent_dir_path_obj

ifelse(!dir.exists(file.path(subdir)),
        dir.create(file.path(subdir),recursive=T),
        "Directory Exists")

figs_dir_path <- file.path(subdir, "figs_postfilterDf")

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

#### =========================================== ####
#### Make BarPlot (post filter-RBC) ####
#### =========================================== ####
DefaultAssay(seurat_obj) <- "RNA"

# Factor sample names
if (is.factor(seurat_obj@meta.data$orig.ident)) {
  seurat_obj@meta.data$orig.ident <- droplevels(seurat_obj@meta.data$orig.ident)
}
if (is.factor(seurat_obj@meta.data$sample)) {
  sample_names <- sort(unique(droplevels(seurat_obj@meta.data$sample)))
} else {
  sample_names <- sort(unique(seurat_obj@meta.data$sample))
}
seurat_obj@meta.data$sample <- factor(
  seurat_obj@meta.data$sample, 
  levels = sample_names
)
sample_names_indexed <- sample_names
names(sample_names_indexed) <- seq(1,length(sample_names))

print("Bar plot Cell Count per Sample")
p <- seurat_obj@meta.data %>%
    ggplot(aes(x=sample, fill=sample)) + 
    labs(x = "", y = "Cell Count") +
    scale_fill_manual(values = sample_palette) + 
    labs(fill="Sample ID") + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=10), legend.position = "none")
ggsave(file.path(figs_dir_path, paste0(filename,"_cellcounts_filtRBC.tiff")), 
      plot = p, units="in", width=size*1, height=size*0.8, dpi=300, compression = 'lzw')

# Table of Cell Counts per Sample
df <- seurat_obj@meta.data %>% 
  group_by(sample) %>% 
  summarise( 
    ncells = n(), 
    noveltyscore.avg = mean(log10GenesPerUMI),
    noveltyscore.sd = sd(log10GenesPerUMI),
    ngene.mad = mad(nGene, constant = 1),
    numi.mad = mad(nUMI, constant = 1),
    mitoRatio.mad = mad(mitoRatio, constant = 1)
  )

# Count number of detected genes in this assay and avg novelty score
# first determine number of genes (rows) not expressed in any cell (column)
obj.list <- SplitObject(seurat_obj, split.by="sample")
obj.list.sorted <- obj.list[order(names(obj.list))]
rm("obj.list")
df$ngenes <- 0

for (i in 1:length(obj.list.sorted)) {
  # Count number of detected genes in this assay and avg novelty score
  # first determine number of genes (rows) not expressed in any cell (column)
  total_num_genes <- nrow(obj.list.sorted[[i]])
  num_undetected_genes <- sum(tabulate(obj.list.sorted[[i]]$RNA@counts@i + 1) == 0) # any zeroes = undetected genes
  df$ngenes[i] <- total_num_genes - num_undetected_genes
  print(paste(
    "There are", num_undetected_genes, "genes out of", dim(obj.list.sorted[[i]])[1], "genes undetected in",  
    names(obj.list.sorted)[i], "resulting in", df$ngenes[i], "detected genes."))
}
rm("obj.list.sorted")
write.csv(df, file= file.path(figs_dir_path, paste0(filename,"_samplecounts_filtRBC.csv")), row.names = F)

## ========================================= ##
## Filter out Doublets ##
## ========================================= ##

DefaultAssay(seurat_obj) <- "RNA"
DefaultAssay(seurat_obj_2) <- "RNA"

meta_1 <- seurat_obj@meta.data['doublet_finder']
meta_1 <- rownames_to_column(meta_1, "row_names")

meta_2 <- seurat_obj_2@meta.data['doublet_finder']
meta_2 <- rownames_to_column(meta_2, "row_names")

print(paste("Number of rows in object 1", nrow(meta_1))) 
print(paste("Number of rows in object 2", nrow(meta_2)))

print("Number of Doublets in obj_1 (before merging labels)")
dplyr::count(meta_1, doublet_finder == "Doublet")
print("Number of Doublets in obj_2 (before merging labels)")
dplyr::count(meta_2, doublet_finder == "Doublet")

## match doublets across GE and GETE objects
meta_1 <- dplyr::left_join (meta_1, meta_2, 
                            by = "row_names",
                            suffix = c(".x", ".y"),
                            unmatched = "error") %>%
          dplyr::mutate(doublets_geXgete = if_else(doublet_finder.x == "Doublet" & 
                                                   doublet_finder.y == "Doublet", 1,         
                                           if_else((doublet_finder.x == "Singlet" & 
                                                   doublet_finder.y == "Doublet") | 
                                                   (doublet_finder.x == "Doublet" & 
                                                   doublet_finder.y == "Singlet"), 
                                                   2, 0) ))

print("Number of Doublets in both objects")
dplyr::count(meta_1, doublets_geXgete ==1)
print("Number of Doublets in one of the objects")
dplyr::count(meta_1, doublets_geXgete ==2)
head(meta_1 %>% dplyr::filter(doublets_geXgete == 2))
print("Number of Singlets in both objects")
dplyr::count(meta_1, doublets_geXgete ==0)

write.csv(meta_1, file.path(figs_dir_path, paste0(filename,"_doublet_summary_geXgete.csv")))

rownames(meta_1) <- meta_1$row_names
meta_1$row_names <- NULL
meta_1$doublet_finder.x <- NULL
meta_1$doublet_finder.y <- NULL
head(meta_1)
seurat_obj <- AddMetaData(seurat_obj, meta_1, 
                          col.name = "doublets_geXgete")

rm("meta_1")
rm("meta_2")
rm("seurat_obj_2")

#### Filter out doublets ####
ncells_old <- dim(seurat_obj)[2]
seurat_obj <- subset(
  x = seurat_obj, 
  subset = (doublets_geXgete == 0)
  # subset = (doublet_finder == "Singlet")
) 
print(paste(ncells_old - dim(seurat_obj)[2], "cells were removed from object"))

#### =========================================== ####
#### Make BarPlot (post filter-Doublets) ####
#### =========================================== ####
DefaultAssay(seurat_obj) <- "RNA"

# Factor sample names
if (is.factor(seurat_obj@meta.data$orig.ident)) {
  seurat_obj@meta.data$orig.ident <- droplevels(seurat_obj@meta.data$orig.ident)
}
if (is.factor(seurat_obj@meta.data$sample)) {
  sample_names <- sort(unique(droplevels(seurat_obj@meta.data$sample)))
} else {
  sample_names <- sort(unique(seurat_obj@meta.data$sample))
}
seurat_obj@meta.data$sample <- factor(
  seurat_obj@meta.data$sample, 
  levels = sample_names
)
sample_names_indexed <- sample_names
names(sample_names_indexed) <- seq(1,length(sample_names))

print("Bar plot Cell Count per Sample")
p <- seurat_obj@meta.data %>%
    ggplot(aes(x=sample, fill=sample)) + 
    labs(x = "", y = "Cell Count") +
    scale_fill_manual(values = sample_palette) + 
    labs(fill="Sample ID") + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=10), legend.position = "none")
ggsave(file.path(figs_dir_path, paste0(filename,"_cellcounts_filtDoublets.tiff")), 
      plot = p, units="in", width=size*1, height=size*0.8, dpi=300, compression = 'lzw')

# Table of Cell Counts per Sample
df <- seurat_obj@meta.data %>% 
  group_by(sample) %>% 
  summarise( 
    ncells = n(), 
    noveltyscore.avg = mean(log10GenesPerUMI),
    noveltyscore.sd = sd(log10GenesPerUMI),
    ngene.mad = mad(nGene, constant = 1),
    numi.mad = mad(nUMI, constant = 1),
    mitoRatio.mad = mad(mitoRatio, constant = 1)
  )

# Count number of detected genes in this assay and avg novelty score
# first determine number of genes (rows) not expressed in any cell (column)
obj.list <- SplitObject(seurat_obj, split.by="sample")
obj.list.sorted <- obj.list[order(names(obj.list))]
rm("obj.list")
df$ngenes <- 0

for (i in 1:length(obj.list.sorted)) {
  # Count number of detected genes in this assay and avg novelty score
  # first determine number of genes (rows) not expressed in any cell (column)
  total_num_genes <- nrow(obj.list.sorted[[i]])
  num_undetected_genes <- sum(tabulate(obj.list.sorted[[i]]$RNA@counts@i + 1) == 0) # any zeroes = undetected genes
  df$ngenes[i] <- total_num_genes - num_undetected_genes
  print(paste(
    "There are", num_undetected_genes, "genes out of", dim(obj.list.sorted[[i]])[1], "genes undetected in",  
    names(obj.list.sorted)[i], "resulting in", df$ngenes[i], "detected genes."))
}
rm("obj.list.sorted")
write.csv(df, file= file.path(figs_dir_path, paste0(filename,"_samplecounts_filtDoublets.csv")), row.names = F)

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
saveRDS(seurat_obj, file = file.path(parent_dir_path,output_dir_name, paste0(filename, "_filtDf_int.rds")) )
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

p <- DimPlot(seurat_obj, reduction = "umap", group.by = "sample") 
ggsave(file.path(figs_dir_path, paste0(filename,"_UMAP_sample_filtDoublets.tiff")), 
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')
print("Exported UMAP by sample")

### Project Minimum number of PCs to UMAP plot
obj_integrated <- RunUMAP(obj_integrated, 
                          reduction = "pca", 
                          dims = 1:min_pc, 
                          umap.method = "uwot", 
                          metric = "cosine")
# saveRDS(obj_integrated, file = file.path(parent_dir_path, output_dir_name, paste0(filename, "_filtDf_umap.rds")) )
# print(paste("ScaleData, RunPCA, and RunUMAP completed in:", parent_dir_path, output_dir_name))
print(paste("RunUMAP completed"))

#### =========================================== ####
#### Cluster: K-nearest neighbor graph
#### =========================================== ####
# rename old cluster annotations
metadata <- seurat_obj@meta.data
metadata <- rename(metadata, preQC_int_snn_res03 = "integrated_snn_res.0.3")
metadata <- rename(metadata, preQC_int_snn_res04 = "integrated_snn_res.0.4")
metadata$seurat_clusters <- NULL
seurat_obj@meta.data <- metadata
glimpse(metadata)
rm("metadata")

## identify clusters
seurat_obj <- FindNeighbors(seurat_obj,dims=1:min_pc,reduction="pca")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
saveRDS(seurat_obj, file = file.path(subdir, paste0(filename, "_filtDf_cluster.rds")))
cat("KNN clustering complete. Saved seurat objects to", subdir,"\n")

#### =========================================== ####
#### Find unique markers per cluster ####
#### =========================================== ####
## Use normalized counts for DGE (Luecken & Theis 2019)
DefaultAssay(seurat_obj) <- "RNA" 
Idents(seurat_obj) <- "integrated_snn_res.0.3"

# Use normalized counts to perform differential gene analysis
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE) # default: LogNormalize
seurat_obj <- ScaleData(seurat_obj,verbose = FALSE)
print("Normalization and Centering + Scaling Complete")

print("Finding Markers...")
markers <- FindAllMarkers(
  seurat_obj, 
  only.pos = FALSE, 
  min.pct = 0.25, 
  logfc.threshold = 0.1)

write.csv(markers, file = file.path(subdir, paste0(filename, "_markers_0.3.csv")), row.names=TRUE)
rm("markers")

## DEG of resolution 0.4 clusters
Idents(seurat_obj) <- "integrated_snn_res.0.4"

## run function
markers <- FindAllMarkers(
  seurat_obj, 
  only.pos = FALSE, 
  min.pct = 0.25, 
  logfc.threshold = 0.1)

write.csv(markers, file = file.path(subdir, paste0(filename, "_markers_0.4.csv")), row.names=TRUE)
rm("markers")

#### End of Script ####
sessionInfo()
