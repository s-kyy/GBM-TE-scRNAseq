
# import 2 files, seurat object + cellular marker CSV (one at a time)

library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least 2 filepaths must be supplied: [xxx.rds] [xxx_markers_0.x.csv]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1]) & file.exists(args[2])) { 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_name_obj <- basename(parent_dir_path_obj)
    marker_path <- args[2]
#     parent_dir_path_marker <- dirname(marker_path)
#     parent_dir_name_marker <- basename(parent_dir_path_marker)
  } else {
    stop("Filepaths provided do not exist. Exiting...", call.=FALSE)
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

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr

set.seed(108)

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
seurat.obj <- readRDS(obj_path) 
size    <- 5
cluster_res03 <- "integrated_snn_res.0.3"
cluster_res04 <- "integrated_snn_res.0.4"

if (grepl("healthy", parent_dir_name_obj, fixed = TRUE)) {
  sample_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#CC79A7", "#ff716e",
    "#999999", "#0072B2", "#194c76", 
    "#D55E00", "#3a4f41", "#6699cc", "#713e5a")
  female_samples <- c("SRR9262922", "SRR9262937",
                        "SRR9264382", "SRR9264383",
                        "SRR9264388")

} else if (grepl("gbm", parent_dir_name_obj, fixed = TRUE)) {
  sample_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#CC79A7", "#ff716e",
    "#999999", "#0072B2", "#194c76")
  female_samples <- "SF11209"
}

figs_dir_path <- file.path(parent_dir_path_obj, "figs")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

print(seurat.obj)
#### =========================================== ####
#### Variable Features Plot ####
#### =========================================== ####
# DefaultAssay(seurat.obj) <- "RNA"
# plot1_ge <- VariableFeaturePlot(
#               object = ge_gbm.sc.intergrated,
#               assay = 
#               )
# plot2_ge <- LabelPoints(plot = plot1_ge, points = top10, repel = TRUE)

#### =========================================== ####
#### UMAP plots by cluster ####
#### =========================================== ####
DefaultAssay(seurat.obj) <-"RNA"

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res03) 
ggsave(file.path(figs_dir_path, paste0(filename, "_", cluster_res03,"_UMAP_nolabels.tiff")),
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res03, 
              label = TRUE, 
              label.size = 2.5) 
ggsave(file.path(figs_dir_path, paste0(filename, "_", cluster_res03,"_UMAP.tiff")),
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res04) 
ggsave(file.path(figs_dir_path, paste0(filename, "_", cluster_res04,"_UMAP_nolabel.tiff")),
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res04, 
              label = TRUE, 
              label.size = 2.5) 
ggsave(file.path(figs_dir_path, paste0(filename, "_", cluster_res04,"_UMAP.tiff")),
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')

print("UMAPs by cluster exported")

#### =========================================== ####
#### UMAP plots by sample ####
#### =========================================== ####

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res03,
              split.by = "sample") 
ggsave(file.path(figs_dir_path, paste0(filename, "_sample_split_", cluster_res03,"_UMAP.tiff")),
       plot = p, units="in", width=size*6, height=size*0.8, dpi=300, compression = 'lzw')

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res04,
              split.by = "sample") 
ggsave(file.path(figs_dir_path, paste0(filename, "_sample_split_", cluster_res04,"_UMAP.tiff")),
       plot = p, units="in", width=size*6, height=size*0.8, dpi=300, compression = 'lzw')

if (grepl("gbm", parent_dir_name_obj, fixed = TRUE)) {
       
  p <- DimPlot(seurat.obj, reduction = "umap", 
       group.by = "sample_orig") 
  ggsave(file.path(figs_dir_path, paste0(filename, "_sample_orig","_UMAP.tiff")),
       plot = p, units="in", width=size*1.3, height=size*1, dpi=300, compression = 'lzw')
  print("UMAPs by sample_orig exported")

  p <- DimPlot(seurat.obj, reduction = "umap", 
                group.by = cluster_res03,
                split.by = "sample_orig") 
  ggsave(file.path(figs_dir_path, paste0(filename, "_sample_orig_split_", cluster_res03,"_UMAP.tiff")),
         plot = p, units="in", width=size*7, height=size*0.8, dpi=300, compression = 'lzw')
  
  p <- DimPlot(seurat.obj, reduction = "umap", 
                group.by = cluster_res04,
                split.by = "sample_orig") 
  ggsave(file.path(figs_dir_path, paste0(filename, "_sample_orig_split_", cluster_res04,"_UMAP.tiff")),
         plot = p, units="in", width=size*7, height=size*0.8, dpi=300, compression = 'lzw')
  
  print("UMAPs split by sample exported")
}

#### =========================================== ####
#### UMAP plots by female sample ####
#### =========================================== ####
# Highlight Female samples
p <- DimPlot(seurat.obj, 
              reduction = "umap", 
              cells.highlight = WhichCells(object = seurat.obj, 
              expression = (sample == female_samples) ))
ggsave(file.path(figs_dir_path, paste0(filename, "_female_samples","_UMAP.tiff")),
      plot = p, units="in", width=size*1.4, height=size*1, dpi=300, compression = 'lzw')

print("UMAPs by female_sample exported")

# Create Heatmap with unbiased clustering of differentially expressed genes 

# Create UMAP plots with Known Markers 

# List of known Markers for each brain cell type

# Filter out genes not expressed in assay

# Create UMAP plots

sessionInfo()