#!usr/bin/env Rscript

.libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least 3 filepaths must be supplied: [gte.rds] [res] [figure_path]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1])) { 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_name_obj <- basename(parent_dir_path_obj)
  } else {
    stop("Filepaths provided do not exist. Exiting...", call.=FALSE)
  }
  
  if (length(args)>=2) {
    # Name of column from integrated dataset 
    cluster_res_num <- as.numeric(args[2])
    cluster_res <- paste0("integrated_snn_res.", cluster_res_num)
    # Custom column name for celltypes in metadata
    cluster_res_symbol <- unlist(strsplit(as.character(cluster_res_num), ".", fixed=TRUE))
    cluster_col <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_celltypes")
    if (grepl("gbm", parent_dir_path_obj, fixed = TRUE)) {
      cluster_col_gsc <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_gsctypes")
      cluster_col_gsc_cc <- paste0(cluster_col_gsc,"_cc")
      cluster_col_gbm_neftel <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_gbmneftel")
      cluster_col_gbm_tcga <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_gbmtcga")
      cluster_col_tumour <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_tumour")
    }
  } else {
    cluster_res <- "integrated_snn_res.0.4"
    cluster_col <- "int04_celltypes"
    if (grepl("gbm", parent_dir_path_obj, fixed = TRUE)) {
      cluster_col_gsc <- "int04_gsctypes"
      cluster_col_gsc_cc <- "int04_gsctypes_cc"
      cluster_col_gbm_neftel <- "int04_gbmneftel"
      cluster_col_gbm_tcga <- "int04_gbmtcga"
      cluster_col_tumour <- "int04_tumour"
    }
  }

  if (length(args)>=3) {
    figs_dir_name <- args[3]
  } else {
    figs_dir_name <- "figs_monocle3"
  }
  

}
#### ===================================================================== ####
#### Import Packages ####
#### ===================================================================== ####

set.seed(108)

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
library(patchwork)
library(magrittr)
library(tidyverse) 
library(dplyr)
library(viridis)

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr
conflict_prefer("rowRanges", "MatrixGenerics") # required in new_cell_data_set()

# dynamic loading of udunits in cedar
udunits_dir = "~/bin/udunits"
dyn.load(file.path(udunits_dir, "local","lib","libudunits2.so.0"))

set.seed(108)

#### ===================================================================== ####
#### Load Datasets ####
#### ===================================================================== ####
seurat_obj <- readRDS(obj_path) 
dim(seurat_obj)
size <- 7

figs_dir_path <- file.path(parent_dir_path_obj, figs_dir_name)

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")


#### ===================================================================== ####
#### Subset out non-tumour and immune cells. ####
#### ===================================================================== ####

seurat_obj <- subset(x = seurat_obj, subset = (!!sym(cluster_col_tumour) == 1) & (!!sym(cluster_col_gsc) != "Tumour Endothelia"))

dim(seurat_obj)
unique(seurat_obj$int06_celltypes)
unique(seurat_obj$int06_gsctypes)
glimpse(seurat_obj@meta.data)
DefaultAssay(seurat_obj) <- "integrated"
Idents(seurat_obj) <- cluster_col
print("Seurat Object subsetted")

#### ===================================================================== ####
#### Perform Monocle Analysis ####
#### ===================================================================== ####

## Extract sparse matrix, metadata and gene labels to create monocle object
cds <- as.cell_data_set(seurat_obj)
cds <- cluster_cells(cds)
print("Monocle object clustered. Exporting UMAPs and saving object...")

# Export cluster and partition figures
p <- plot_cells(cds, show_trajectory_graph = FALSE, group_label_size = 4) + 
  theme(axis.line=element_blank(),
    axis.text.x=element_blank(),axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),axis.title.y=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
ggsave(file.path(figs_dir_path, paste0(filename, "_monocleUMAPS_clusters.tiff")), 
  plot = p, units="in", width=size, height=size, dpi=300, compression = 'lzw')

p <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE, group_label_size = 4) +
  theme(axis.line=element_blank(),
    axis.text.x=element_blank(),axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),axis.title.y=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())
ggsave(file.path(figs_dir_path, paste0(filename, "_monocleUMAPS_partitions.tiff")), 
  plot = p, units="in", width=size, height=size, dpi=300, compression = 'lzw')

saveRDS(cds, file=file.path(parent_dir_path_obj, paste0(filename, "_monocle.rds")))

### End of Script
sessionInfo()
