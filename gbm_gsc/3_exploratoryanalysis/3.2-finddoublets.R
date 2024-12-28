#!usr/bin/env Rscript

.libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 1 filepath must be supplied: [xxx.rds] [sample_name]", call.=FALSE)
} else if (length(args)<2) {
  # verify filepaths
  if (file.exists(args[1])){ 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    sample_name <- ""    
    parent_dir_path_obj <- dirname(dirname(obj_path))
    parent_dir_name_obj <- basename(parent_dir_path_obj)
  } else {
    stop("Filepath provided does not exist. Exiting...", call.=FALSE)
  }
} else if (length(args)<=2) {
  if (file.exists(args[1])){ 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    sample_name <- args[2]
    parent_dir_path_obj <- dirname(dirname(obj_path))
    parent_dir_name_obj <- basename(parent_dir_path_obj)
  } else {
    stop("Filepath provided does not exist. Exiting...", call.=FALSE)
  }
}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 
library(dplyr)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)

# DoubletFinder@03e9f37f891ef76a23cc55ea69f940c536ae8f9f (April 10, 2024)
# remotes::install_github(repo='chris-mcginnis-ucsf/DoubletFinder', ref = remotes::github_pull(176))
library(DoubletFinder)

set.seed(108)
options(warn=1) #print warning messages as they occur

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
# Ensure objects have reached the RunUMAP step
seurat.obj <- readRDS(obj_path) 
size    <- 5
ndims <- 25

# multiplet_rate_10x <- data.frame(
#   multiplet_rates = c(0.008,0.016,0.024,0.032,0.04,0.048,0.056,
#                         0.064,0.072,0.08,0.096,0.112,0.128,
#                         0.144,0.16,0.176,0.192,0.208,0.224,0.24),
#   loaded_cells =    c(3130,6320,9550,12800,16100,19500,
#                         22900,26300,29800,33300,40500,47800,55400,
#                         63200,71200,79500,88100,96900,10600,11500),
#   recovered_cells = c(2000,4000,6000,8000,10000,12000,14000,
#                         16000,18000,20000,24000,28000,32000,36000,
#                         40000,44000,48000,52000,56000,60000)
# )

print("Objects loaded")

subdir <- file.path(parent_dir_path_obj, paste0(format(Sys.Date(), "%Y%m%d"), "_", sample_name,"_", filename, "_doublets"))

ifelse(!dir.exists(file.path(subdir)),
        dir.create(file.path(subdir),recursive=T),
        "Directory Exists")

figs_dir_path <- file.path(subdir, "figs_doublets")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

## ========================================= ##
## Function: Preprocess seurat subset ##
## ========================================= ##
preprocess_seurat_subset <- function(obj, ndims) {
  
  ## STEPS to re-cluster after subsetting
  # Find variable features in the subset using the RNA assay
  # Run ScaleData on the integrate assay on the new set of variable features
  # Run PCA on the integrated assay using the new set of variable features
  # Run FindNeighbors and FindClusters using the new PC dimensions
  # source: https://github.com/satijalab/seurat/issues/5532

  #### FindVariableFeatures + Scale Data ~5min, <200Mb ####
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = TRUE)
  obj <- FindVariableFeatures(obj, 
                              selection.method="vst",
                              nfeatures=2500, #default=2000
                              verbose = TRUE)

  #### Create Variable Feature Plots ####
  top10 <- head(VariableFeatures(obj), 10)
  p <- VariableFeaturePlot(obj)
  p <- LabelPoints(plot = p, points = top10, repel = TRUE)
  ggsave(file.path(figs_dir_path, paste0(filename,"_varfeatgenes.tiff")), 
        plot = p, units="in", width=size*1.1, height=size*0.8, dpi=300, compression = 'lzw')

  # DefaultAssay(obj) <- "integrated"
  obj <- ScaleData(object = obj, verbose = FALSE, features = VariableFeatures(obj[["RNA"]])) 

  # Run PCA
  obj <- RunPCA(obj, 
                npcs = ndims, 
                verbose = FALSE)
  # saveRDS(obj, file = file.path(subdir, paste0(filename, "_filt.rds")) )
  # print(paste("ScaleData, RunPCA completed in:", subdir))

  #### Test PCA levels ####
  # Determine percent of variation associated with each PC
  # DefaultAssay(obj) <- "integrated"
  pct_var_per_pc <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100

  # Calculate cumulative percents for each PC
  cum_pct_per_pc <- cumsum(pct_var_per_pc)

  # Determine which PC exhibits a cumulative percentage of variation 
  # greater than 90% and variation associated with the PC is less than 5%
  pc_most_var <- which(cum_pct_per_pc > 90 & pct_var_per_pc < 5)[1]
  print(paste("Minimum PC that retains more than 90% variation and less than 5% variation compared to the next PC:", pc_most_var))

  # Determine the difference between variation of PC and subsequent PC
  pc_10 <- sort(which((pct_var_per_pc[1:length(pct_var_per_pc) - 1] - pct_var_per_pc[2:length(pct_var_per_pc)]) > 0.1), decreasing = T)[1] + 1

  # last point where change of % of variation is more than 0.1%.
  print(paste("Minimum PC with a difference in variation of 10% compared to next PC:", pc_10))
  co2

  # Minimum of the two calculation
  min_pc <- min(pc_most_var, pc_10)
  print(paste("Minumum PC between the two options:", min_pc))

  plot_df <- data.frame(dimensions = 1:length(pct_var_per_pc),
            stdev = obj[["pca"]]@stdev,
            pct_var_per_pc = pct_var_per_pc,
            cum_pct_per_pc = cum_pct_per_pc)
  print(plot_df)
  write.csv(plot_df, file.path(figs_dir_path, paste0(filename, "_pca.csv") ))

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
  ggsave(file.path(figs_dir_path, paste0(filename,"_elbow.tiff")), 
    plot = p, units="in", width=size*0.7, height=size*0.7, dpi=300, compression = 'lzw')
  print("Exported ElbowPlot")

  obj <- RunUMAP(obj, 
                  reduction = "pca", 
                  dims = 1:min_pc, 
                  umap.method = "uwot", 
                  metric = "cosine")
  saveRDS(obj, file = file.path(subdir, paste0(filename, "_filt_umap_temp.rds")) )
  print(paste("RunUMAP completed in:", subdir))

  #### Output UMAP plots by sample ####
  ggsave(file.path(figs_dir_path, paste0(filename,"_UMAP_sample.tiff")), 
        plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')
  print("Exported UMAP by sample")

  #### Cluster: K-nearest neighbor graph
  obj <- FindNeighbors(obj,dims=1:min_pc,reduction="pca")
  obj <- FindClusters(obj, resolution = 0.3)
  obj <- FindClusters(obj, resolution = 0.4)
  colnames(obj@meta.data)
  # saveRDS(obj, file = file.path(subdir, paste0(filename, "_clustered.rds")))
  print("KNN clustering complete.")

  meta <- obj@meta.data

  meta$`integrated_snn_res.0.3` <- as.character(meta$`integrated_snn_res.0.3`)
  meta$`integrated_snn_res.0.4` <- as.character(meta$`integrated_snn_res.0.4`)
  meta$`RNA_snn_res.0.3` <- as.character(meta$`RNA_snn_res.0.3`)
  meta$`RNA_snn_res.0.4` <- as.character(meta$`RNA_snn_res.0.4`)
  meta$seurat_clusters <- as.character(meta$seurat_clusters)

  obj@meta.data <- meta

  results <- list("obj" = obj, "min_pc" = min_pc)
  return(results)
}

#### =========================================== ####
#### Function: findDoublets ####
#### =========================================== ####
#### https://github.com/chris-mcginnis-ucsf/DoubletFinder

findDoublets <- function(seu, dims, rate) {

  print(paste("Preprocess Seurat Object"))

  results <- preprocess_seurat_subset(seu, dims)
  obj_sample <- results[["obj"]]
  dims <- results[["min_pc"]]
  rm("results")

  dplyr::glimpse(obj_sample@meta.data)

  #### pK Identification (no ground truth) - number of expected nearest neighbors that are multiplets (neighborhood size) 
  sweep.res.list <- paramSweep(obj_sample, PCs = 1:dims, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn_sample <- find.pK(sweep.stats)
  write.csv(bcmvn_sample, file=file.path(figs_dir_path, paste0(sample_name, "_",filename,"_bcmvn.csv")))

  #### Take the max value in bimodality coefficient (BCmvn) distribution. 
  optimal_pK <- bcmvn_sample %>%
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal_pK <- as.numeric(as.character(optimal_pK[[1]]))
  print(paste("Optimal pK threshold for sample", unique(obj_sample@meta.data$sample), "is"))
  print(optimal_pK)

  #### nEXP - number of expected real multiplets 
  cluster_annotations <- as.character(obj_sample@meta.data$`RNA_snn_res.0.4`)
  if (!is.null(cluster_annotations) && length(cluster_annotations) > 0) {
    homotypic_prop <- modelHomotypic(cluster_annotations)
    print(table(cluster_annotations))
    print(paste("Homotypic Proportion for sample", unique(obj_sample@meta.data$sample)))
    print(homotypic_prop)
    nEXP_poi <- round(rate * nrow(obj_sample@meta.data))
    nEXP_poi_adj <- round(nEXP_poi * (1 - homotypic_prop))
    print(paste("nEXP_poi_adj (expected real multiplets) =", nEXP_poi_adj))

    #### Run DoubletFinder 
    obj_sample <- doubletFinder(seu = obj_sample, PCs = 1:dims, 
                            pK = optimal_pK, nExp = nEXP_poi_adj) #pN = 0.25 (default)

    #### rename doublet finder annotations column and extract metadata 
    meta <- obj_sample@meta.data

    colnames(meta)[grepl('DF.classifications.*', colnames(meta))] <- "doublet_finder"
    results <- meta['doublet_finder']
    results <- rownames_to_column(results,"row_names")
    print(head(results))
    return(results)

  } else {
    stop("No Cluster annotations available. Exiting...")
  }
}

#### =========================================== ####
#### Main ####
#### =========================================== ####

#### Filter RBC cells ####
DefaultAssay(seurat.obj) <- "RNA"
filtered_obj <- subset(
  x=seurat.obj, 
  subset= 
  (HBA1 < 1) &
  (HBA2 < 1) &
  (HBB < 1) &
  (HBG2 < 1) &
  (HBZ < 1) & 
  (HBM < 1) & 
  (HBD < 1) & 
  (HBE1 < 1) &
  (HBG1 < 1) &
  (HBQ1 < 1)  
) 
print(paste(dim(seurat.obj)[2] - dim(filtered_obj)[2], "cells were removed from object"))
rm("seurat.obj")

#### Split objects ####
if (sample_name == "healthy") {

  # split healthy dataset by sample. 
  dplyr::glimpse(filtered_obj@meta.data)
  filtered_obj@meta.data$sample <- as.character(filtered_obj@meta.data$sample)
  obj_list <- SplitObject(object = filtered_obj, split.by="sample")
  multiplet_rate_10x <- 0.01
  
  #### Run DoubletFinder ####
  doubletfinder_results <- vector(mode='list', length=length(obj_list))
  
  for (i in 1:length(obj_list)) { 
    doubletfinder_results[[i]] <- findDoublets(seu = obj_list[[i]], 
                                                dims = ndims, 
                                                rate = multiplet_rate_10x)
  }

} else if (sample_name == "gbm") {
  
  # copy GBM_WANG2020 sample names to sample_orig column
  dplyr::glimpse(filtered_obj@meta.data)
  filtered_obj@meta.data$sample_orig <- as.character(filtered_obj@meta.data$sample_orig)
  filtered_obj@meta.data$sample <- as.character(filtered_obj@meta.data$sample)

  for (i in which(is.na(filtered_obj@meta.data$sample_orig))){
    filtered_obj@meta.data$sample_orig[i] <- as.character(filtered_obj@meta.data$sample[i])
  }
  p <- DimPlot(filtered_obj, reduction = "umap", group.by = "sample_orig") 
  ggsave(file.path(figs_dir_path, paste0(filename,"_UMAP_sample_orig.tiff")), 
        plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')
  print("Exported UMAP by sample_orig")

  # split gbm_bhaduri sample by lane and gbm_wang by sample
  obj_list <- SplitObject(object = filtered_obj, split.by="sample_orig")
  multiplet_rate_10x <- rep(c(0.008, 0.026), c(12,3))
  
  #### Run DoubletFinder ####
  doubletfinder_results <- vector(mode='list', length=length(obj_list))

  for (i in 1:length(obj_list)) { 
    doubletfinder_results[[i]] <- findDoublets(seu = obj_list[[i]], 
                                                dims = ndims, 
                                                rate = multiplet_rate_10x[i])
  }
}
print("Doublets annotated")

#### =========================================== ####
#### Merge Metadata ####
#### =========================================== ####
# merge dataframes from doubletfinder output
doubletfinder_results_all <- data.frame(bind_rows(doubletfinder_results))
rownames(doubletfinder_results_all) <- doubletfinder_results_all$row_names
doubletfinder_results_all$row_names <- NULL
head(doubletfinder_results_all)
filtered_obj <- AddMetaData(filtered_obj, doubletfinder_results_all, 
                          col.name = 'doublet_finder')
print("Doublet Finder results metadata added to seurat object")

#### =========================================== ####
#### Export Object
#### =========================================== ####
saveRDS(filtered_obj, file = file.path(subdir, paste0(filename, "_ANNdoublets.rds")) )
print(paste("Seurat Object exported to",subdir))

#### =========================================== ####
#### Figures
#### =========================================== ####
# Reference: https://youtu.be/lwOGXR7kRZ0?si=o-6WqVuz1l1g25X2

if (sample_name == "healthy") {

  print(colnames(filtered_obj@meta.data))

  p <- ggplot(filtered_obj@meta.data, 
              aes (x = sample, y = nGene, fill = doublet_finder, 
                group = interaction(sample, doublet_finder))) + 
              geom_violin(width = 0.8, position = position_dodge(width = 0.8), alpha = 0.2, scale = "width") +
              labs(x = "", y = "Number of Genes Detected") +
              theme(
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12, color="black", ),
                axis.title.y.left = element_text(size=12, face="bold"),
                axis.title.x=element_blank(),
                axis.line = element_line(linewidth=0.5, colour="black"),
                panel.background=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank(),
                legend.position = 'right',
                legend.title = element_blank()) 
  p <- p + geom_point(
            mapping = aes(color = doublet_finder, group = interaction(sample, doublet_finder)), alpha=0.9, size=0.2, position = position_jitterdodge(jitter.height=0.25, jitter.width=0.25,dodge.width=0.8)) + 
          geom_boxplot(fill="white", outlier.colour=NA, width = 0.2, position=position_dodge(width=0.8))
  ggsave(file.path(figs_dir_path, paste0("Vlnplt_nGene_doubletfinder_sample.tiff")), 
        plot = p, units="in", width=size*1.3, height=size*1.2, dpi=300, compression = 'lzw')

  p <- ggplot(filtered_obj@meta.data, 
              aes (x = sample, y = nUMI, fill = doublet_finder, 
                group = interaction(sample, doublet_finder))) + 
              geom_violin(width = 0.8, position = position_dodge(width = 0.8), alpha = 0.2, scale = "width") +
              labs(x = "", y = "Number of Genes Detected") +
              theme(
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12, color="black", ),
                axis.title.y.left = element_text(size=12, face="bold"),
                axis.title.x=element_blank(),
                axis.line = element_line(linewidth=0.5, colour="black"),
                panel.background=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank(),
                legend.position = 'right',
                legend.title = element_blank()) 
  p <- p + geom_point(
            mapping = aes(color = doublet_finder, group = interaction(sample, doublet_finder)), alpha=0.9, size=0.2, position = position_jitterdodge(jitter.height=0.25, jitter.width=0.25,dodge.width=0.8)) + 
          geom_boxplot(fill="white", outlier.colour=NA, width = 0.2, position=position_dodge(width=0.8))
  ggsave(file.path(figs_dir_path, paste0("Vlnplt_nUMI_doubletfinder_sample.tiff")), 
        plot = p, units="in", width=size*1.3, height=size*1.2, dpi=300, compression = 'lzw')

  p <- ggplot(filtered_obj@meta.data, 
              aes (x = sample, y = mitoRatio, fill = doublet_finder, 
                group = interaction(sample, doublet_finder))) + 
              geom_violin(width = 0.8, position = position_dodge(width = 0.8), alpha = 0.2, scale = "width") +
              labs(x = "", y = "Number of Genes Detected") +
              theme(
                axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12, color="black", ),
                axis.title.y.left = element_text(size=12, face="bold"),
                axis.title.x=element_blank(),
                axis.line = element_line(linewidth=0.5, colour="black"),
                panel.background=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank(),
                legend.position = 'right',
                legend.title = element_blank()) 
p <- p + geom_point(
          mapping = aes(y = mitoRatio, color = doublet_finder, group = interaction(sample, doublet_finder)), 
          alpha=0.9, size=0.2, position = position_jitterdodge(jitter.height=0, jitter.width=0.05,dodge.width=0.8)) + 
         geom_boxplot(fill="white", outlier.colour=NA, width = 0.2, position=position_dodge(width=0.8))
  ggsave(file.path(figs_dir_path, paste0("Vlnplt_mitoRatio_doubletfinder_sample.tiff")), 
        plot = p, units="in", width=size*1.3, height=size*1.2, dpi=300, compression = 'lzw')
  print("Doublets annotated in ViolinPlot")

  #### Summary table
  doublets_summary <- filtered_obj@meta.data %>%
    group_by(sample, doublet_finder) %>%
    summarize(total_count = n(), .groups = 'drop') %>%
    as.data.frame() %>% ungroup() %>%
    group_by(sample) %>%
    mutate(countT = sum(total_count)) %>%
    group_by(doublet_finder, .add = TRUE) %>%
    mutate('percent (%)' = round(100*total_count/countT, 2)) %>%
    dplyr::select(-countT)
  head(doublets_summary)
  write.csv(doublets_summary, file = file.path(figs_dir_path, "doubletfinder_summary.csv"))
  print("csv exported")

} else if (sample_name == "gbm") {

  print(colnames(filtered_obj@meta.data))

  p <- ggplot(filtered_obj@meta.data, 
                  aes (x = sample_orig, y = nGene, fill = doublet_finder, 
                    group = interaction(sample_orig, doublet_finder))) + 
                  geom_violin(width = 0.8, alpha = 0.2, scale = "width",
                    position = position_dodge(width = 0.8)) +
                  labs(x = "", y = "Number of Genes Detected") +
                  theme(
                    axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12, color="black", ),
                    axis.title.y.left = element_text(size=12, face="bold"),
                    axis.title.x=element_blank(),
                    axis.line = element_line(linewidth=0.5, colour="black"),
                    panel.background=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    plot.background=element_blank(),
                    legend.position = 'right',
                    legend.title = element_blank()) 
  p <- p + geom_point(
            mapping = aes(color = doublet_finder, group = interaction(sample_orig, doublet_finder)), 
            alpha=0.9, size=0.2, position = position_jitterdodge(jitter.height=0.25, jitter.width=0.25,dodge.width=0.8)) + 
          geom_boxplot(fill="white", outlier.colour=NA, width = 0.2, position=position_dodge(width=0.8))
  ggsave(file.path(figs_dir_path, paste0("Vlnplt_nGene_doubletfinder_sample_orig.tiff")), 
        plot = p, units="in", width=size*1.3, height=size*1.2, dpi=300, compression = 'lzw')

  p <- ggplot(filtered_obj@meta.data, 
                  aes (x = sample_orig, y = nUMI, fill = doublet_finder, 
                    group = interaction(sample_orig, doublet_finder))) + 
                  geom_violin(width = 0.8, alpha = 0.2, scale = "width",
                    position = position_dodge(width = 0.8)) +
                  labs(x = "", y = "Number of Genes Detected") +
                  theme(
                    axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12, color="black", ),
                    axis.title.y.left = element_text(size=12, face="bold"),
                    axis.title.x=element_blank(),
                    axis.line = element_line(linewidth=0.5, colour="black"),
                    panel.background=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    plot.background=element_blank(),
                    legend.position = 'right',
                    legend.title = element_blank()) 
  p <- p + geom_point(
            mapping = aes(color = doublet_finder, group = interaction(sample_orig, doublet_finder)), 
            alpha=0.9, size=0.2, position = position_jitterdodge(jitter.height=0.25, jitter.width=0.25,dodge.width=0.8)) + 
          geom_boxplot(fill="white", outlier.colour=NA, width = 0.2, position=position_dodge(width=0.8))
  ggsave(file.path(figs_dir_path, paste0("Vlnplt_nUMI_doubletfinder_sample_orig.tiff")), 
        plot = p, units="in", width=size*1.3, height=size*1.2, dpi=300, compression = 'lzw')

  p <- ggplot(filtered_obj@meta.data, 
                  aes (x = sample_orig, y = mitoRatio, fill = doublet_finder, 
                    group = interaction(sample_orig, doublet_finder))) + 
                  geom_violin(width = 0.8, alpha = 0.2, scale = "width",
                    position = position_dodge(width = 0.8)) +
                  labs(x = "", y = "Number of Genes Detected") +
                  theme(
                    axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=12, color="black", ),
                    axis.title.y.left = element_text(size=12, face="bold"),
                    axis.title.x=element_blank(),
                    axis.line = element_line(linewidth=0.5, colour="black"),
                    panel.background=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    plot.background=element_blank(),
                    legend.position = 'right',
                    legend.title = element_blank()) 
  p <- p + geom_point(
            mapping = aes(y = mitoRatio, color = doublet_finder, group = interaction(sample_orig, doublet_finder)), 
            alpha=0.9, size=0.2, position = position_jitterdodge(jitter.height=0, jitter.width=0.05,dodge.width=0.8)) + 
          geom_boxplot(fill="white", outlier.colour=NA, width = 0.2, position=position_dodge(width=0.8))
  ggsave(file.path(figs_dir_path, paste0("Vlnplt_mitoRatio_doubletfinder_sample_orig.tiff")), 
        plot = p, units="in", width=size*1.3, height=size*1.2, dpi=300, compression = 'lzw')  
  print("Doublets annotated in UMAP")

  #### Summary table
  doublets_summary <- filtered_obj@meta.data %>%
    group_by(sample_orig, doublet_finder) %>%
    summarize(total_count = n(), .groups = 'drop') %>%
    as.data.frame() %>% ungroup() %>%
    group_by(sample_orig) %>%
    mutate(countT = sum(total_count)) %>%
    group_by(doublet_finder, .add = TRUE) %>%
    mutate(percent = paste0(round(100*total_count/countT, 2), '%')) %>%
    dplyr::select(-countT)
  head(doublets_summary)
  write.table(doublets_summary, file = file.path(figs_dir_path, paste0(filename, "_doubletfinder_summary.csv")))
  print("DoubletFinder summary csv exported")

}

p <- DimPlot(filtered_obj, reduction = "umap", label =F, 
              cells.highlight = WhichCells(object = filtered_obj, 
                expression = (doublet_finder == "Doublet") )) + NoLegend()
p <- p + theme(axis.line=element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave(file.path(figs_dir_path, paste0("UMAP_doublets.tiff")),
    plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')

#### End of Script ####
sessionInfo()