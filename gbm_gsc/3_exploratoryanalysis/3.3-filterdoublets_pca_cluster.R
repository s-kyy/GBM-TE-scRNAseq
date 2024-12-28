#!usr/bin/env Rscript

#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 1 filepath must be supplied: [xxx.rds] [sample_name]", call.=FALSE)
} else if (length(args)>=1) {
  
  # verify filepaths
  if (file.exists(args[1])){ 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(dirname(obj_path))
    parent_dir_name_obj <- basename(parent_dir_path_obj)
  } else {
    stop("Filepath provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>1) {
    sample_name <- args[2]
  } else {
    sample_name <- ""
  }

} else {
    stop("Error. ", call.=FALSE)
}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/tcga-gbm-R4/renv/library/R-4.0/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()

set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 

set.seed(108)

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
seurat.obj <- readRDS(obj_path) 
size    <- 5
ndims <- 25

subdir <- file.path(getwd(), paste0(format(Sys.Date(), "%Y%m%d"), "_", sample_name,"_", filename))

ifelse(!dir.exists(file.path(subdir)),
        dir.create(file.path(subdir),recursive=T),
        "Directory Exists")

figs_dir_path <- file.path(subdir, "figs_post_filterRBC_postfilterDoublets_pca")

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
DefaultAssay(seurat.obj) <- "RNA"

# Factor sample names
if (is.factor(seurat.obj@meta.data$orig.ident)) {
  seurat.obj@meta.data$orig.ident <- droplevels(seurat.obj@meta.data$orig.ident)
}
if (is.factor(seurat.obj@meta.data$sample)) {
  sample_names <- sort(unique(droplevels(seurat.obj@meta.data$sample)))
} else {
  sample_names <- sort(unique(seurat.obj@meta.data$sample))
}
seurat.obj@meta.data$sample <- factor(
  seurat.obj@meta.data$sample, 
  levels = sample_names
)
sample_names_indexed <- sample_names
names(sample_names_indexed) <- seq(1,length(sample_names))

print("Bar plot Cell Count per Sample")
p <- seurat.obj@meta.data %>%
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
df <- seurat.obj@meta.data %>% 
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
obj.list <- SplitObject(seurat.obj, split.by="sample")
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

#### Filter RBC cells ####
DefaultAssay(seurat.obj) <- "RNA"
filtered_obj <- subset(
  x = seurat.obj, 
  subset = (doublet_finder == "Singlet")
) 
print(paste(dim(seurat.obj)[2] - dim(filtered_obj)[2], "cells were removed from object"))
seurat.obj <- filtered_obj
rm("filtered_obj")

#### =========================================== ####
#### Make BarPlot (post filter-Doublets) ####
#### =========================================== ####
DefaultAssay(seurat.obj) <- "RNA"

# Factor sample names
if (is.factor(seurat.obj@meta.data$orig.ident)) {
  seurat.obj@meta.data$orig.ident <- droplevels(seurat.obj@meta.data$orig.ident)
}
if (is.factor(seurat.obj@meta.data$sample)) {
  sample_names <- sort(unique(droplevels(seurat.obj@meta.data$sample)))
} else {
  sample_names <- sort(unique(seurat.obj@meta.data$sample))
}
seurat.obj@meta.data$sample <- factor(
  seurat.obj@meta.data$sample, 
  levels = sample_names
)
sample_names_indexed <- sample_names
names(sample_names_indexed) <- seq(1,length(sample_names))

print("Bar plot Cell Count per Sample")
p <- seurat.obj@meta.data %>%
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
df <- seurat.obj@meta.data %>% 
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
obj.list <- SplitObject(seurat.obj, split.by="sample")
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
## Scale and Reduce Dimensionality without regression ##
## ========================================= ##
DefaultAssay(seurat.obj) <- "integrated"
# Scale Data ~5min, <200Mb
seurat.obj <- ScaleData(object = seurat.obj, verbose = FALSE) 

seurat.obj@meta.data$sample <- factor(seurat.obj@meta.data$sample, 
                                        levels = sort(unique(seurat.obj@meta.data$sample)))
# Run PCA and UMAP ~20min
seurat.obj <- RunPCA(seurat.obj, 
                      npcs = ndims, 
                      verbose = FALSE)

#### =========================================== ####
#### Test PCA levels ####
#### =========================================== ####
# Determine percent of variation associated with each PC
DefaultAssay(seurat.obj) <- "integrated"
pct_var_per_pc <- seurat.obj[["pca"]]@stdev / sum(seurat.obj[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cum_pct_per_pc <- cumsum(pct_var_per_pc)

# Determine which PC exhibits a cumulative percentage of variation 
# greater than 90% and variation associated with the PC is less than 5%
pc_most_var <- which(cum_pct_per_pc > 90 & pct_var_per_pc < 5)[1]
print(paste("Minimum PC that retains more than 90% variation and less than 5% variation compared to the next PC:", pc_most_var))

# Determine the difference between variation of PC and subsequent PC
pc_10 <- sort(which((pct_var_per_pc[1:length(pct_var_per_pc) - 1] - pct_var_per_pc[2:length(pct_var_per_pc)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
print(paste("Minimum PC with a difference in variation of 0.1% compared to next PC:", pc_10))

# Minimum of the two calculation
min_pc <- min(pc_most_var, pc_10)
print(paste("Minumum PC between the two options:", min_pc))

plot_df <- data.frame(dimensions = 1:length(pct_var_per_pc),
           stdev = seurat.obj[["pca"]]@stdev,
           pct_var_per_pc = pct_var_per_pc,
           cum_pct_per_pc = cum_pct_per_pc)
print(plot_df)
write.csv(plot_df, file.path(figs_dir_path, paste0(filename, "_pca_filtDoublets.csv") ))

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
ggsave(file.path(figs_dir_path, paste0(filename,"_elbow_filtDoublets.tiff")), 
  plot = p, units="in", width=size*0.7, height=size*0.7, dpi=300, compression = 'lzw')
print("Exported ElbowPlot")

#### =========================================== ####
#### Output UMAP plots by sample ####
#### =========================================== ####

p <- DimPlot(seurat.obj, reduction = "umap", group.by = "sample") 
ggsave(file.path(figs_dir_path, paste0(filename,"_UMAP_sample_filtDoublets.tiff")), 
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')
print("Exported UMAP by sample")

#### =========================================== ####
#### Cluster: K-nearest neighbor graph
#### =========================================== ####
seurat.obj <- FindNeighbors(seurat.obj,dims=1:min_pc,reduction="pca")
seurat.obj <- FindClusters(seurat.obj, resolution = 0.3)
seurat.obj <- FindClusters(seurat.obj, resolution = 0.4)
saveRDS(seurat.obj, file = file.path(subdir, paste0(filename, "_clustered_filtDoublets.rds")))
cat("KNN clustering complete. Saved seurat objects to", subdir,"\n")

#### =========================================== ####
#### Find unique markers per cluster ####
#### =========================================== ####
## Use normalized counts for DGE (Luecken & Theis 2019)
DefaultAssay(seurat.obj) <- "RNA" 
Idents(seurat.obj) <- "integrated_snn_res.0.3"

# Use normalized counts to perform differential gene analysis
seurat.obj <- NormalizeData(seurat.obj, verbose = TRUE) # default: LogNormalize

markers <- FindAllMarkers(
  seurat.obj, 
  only.pos = FALSE, 
  min.pct = 0.25, 
  logfc.threshold = 0.1)

write.csv(markers, file = file.path(subdir, paste0(filename, "_markers_0.3.csv")), row.names=TRUE)
rm("markers")

## DEG of resolution 0.4 clusters
Idents(seurat.obj) <- "integrated_snn_res.0.4"

## run function
markers <- FindAllMarkers(
  seurat.obj, 
  only.pos = FALSE, 
  min.pct = 0.25, 
  logfc.threshold = 0.1)

write.csv(markers, file = file.path(subdir, paste0(filename, "_markers_0.4.csv")), row.names=TRUE)
rm("markers")

#### End of Script ####
sessionInfo()
