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
cds <- readRDS(obj_path) 
dim(cds)
size <- 7

figs_dir_path <- file.path(parent_dir_path_obj, figs_dir_name)

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")


#### ===================================================================== ####
#### Monocle Figures ####
#### ===================================================================== ####

# Select Island to study the trajectory on:
partitioned_obj <- subset(as.Seurat(cds), monocle3_partitions == 1)

# Convert back to monocle object and learn the graph
partitioned_obj <- as.cell_data_set(partitioned_obj)
partitioned_obj <- learn_graph(partitioned_obj)
saveRDS(partitioned_obj, file=file.path(parent_dir_path_obj, paste0(filename, "_graph.rds")))
rm("cds")

# Export trajectory graph (by celltype)
p <- partitioned_obj %>% plot_cells(
  color_cells_by = cluster_col, # colour by celltype
  label_cell_groups=FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  # graph_label_size = 3,
  # trajectory_graph_color = "white"
)
ggsave(file.path(figs_dir_path, paste0(filename, "_graph_",cluster_col,".tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

# Export trajectory graph (by celltype) with labels
p <- partitioned_obj %>% plot_cells(
  color_cells_by = cluster_col, # colour by celltype
  label_cell_groups=FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  graph_label_size = 3
)
ggsave(file.path(figs_dir_path, paste0(filename, "_graph_",cluster_col,"_labels.tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

# Export trajectory graph (by gsctypes)
p <- partitioned_obj %>% plot_cells(
  color_cells_by = cluster_col_gsc, # colour by celltype
  label_cell_groups=FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
  # graph_label_size = 3,
  # trajectory_graph_color = "white"
)
ggsave(file.path(figs_dir_path, paste0(filename, "_graph_",cluster_col_gsc,".tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

# Export trajectory graph (by celltype) with labels
p <- partitioned_obj %>% plot_cells(
  color_cells_by = cluster_col_gsc, # colour by celltype
  label_cell_groups=FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  graph_label_size = 3
)
ggsave(file.path(figs_dir_path, paste0(filename, "_graph_",cluster_col_gsc,"_labels.tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

# Export trajectory graph (by gbm_neftel)
p <- partitioned_obj %>% plot_cells(
  color_cells_by = cluster_col_gbm_neftel, # colour by celltype
  label_cell_groups=FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
  # graph_label_size = 3,
  # trajectory_graph_color = "white"
)
ggsave(file.path(figs_dir_path, paste0(filename, "_graph_",cluster_col_gbm_neftel,".tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

# Export trajectory graph (by celltype) with labels
p <- partitioned_obj %>% plot_cells(
  color_cells_by = cluster_col_gbm_neftel, # colour by celltype
  label_cell_groups=FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  graph_label_size = 3
)
ggsave(file.path(figs_dir_path, paste0(filename, "_graph_",cluster_col_gbm_neftel,"_labels.tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

#### ===================================================================== ####
#### Compute Pseudotime: Select root cells by cluster ####
#### ===================================================================== ####
seurat_obj <- as.Seurat(partitioned_obj)
root_cells <- head(rownames(seurat_obj@meta.data[which(FetchData(seurat_obj,cluster_col) == "Mixed Cycling"),]))
rm("seurat_obj")

## Order Cells
partitioned_obj_ordered <- order_cells(partitioned_obj, root_cells = root_cells)

## Export pseudotime graph
p <- partitioned_obj_ordered %>% plot_cells(
  color_cells_by = "pseudotime", 
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE,
  label_roots = FALSE,
  trajectory_graph_color = "black"
  # graph_label_size = 5,
)
ggsave(file.path(figs_dir_path, paste0(filename, "_pseudotime.tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

## Export pseudotime graph (labeled)
p <- partitioned_obj_ordered %>% plot_cells(
  color_cells_by = "pseudotime", 
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE,
  label_roots = TRUE,
  trajectory_graph_color = "black",
  graph_label_size = 3
)
ggsave(file.path(figs_dir_path, paste0(filename, "_pseudotime_labelroot.tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

## Export pseudotime graph (labeled)
p <- partitioned_obj_ordered %>% plot_cells(
  color_cells_by = "pseudotime", 
  label_cell_groups = FALSE, 
  label_leaves = TRUE, 
  label_branch_points = FALSE,
  label_roots = FALSE,
  trajectory_graph_color = "black",
  graph_label_size = 3
)
ggsave(file.path(figs_dir_path, paste0(filename, "_pseudotime_labelleaves.tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

# save pseudotime ordering
pseudotime <- pseudotime(partitioned_obj_ordered) 	
partitioned_obj_ordered@colData$pseudotime <- pseudotime 	

glimpse(partitioned_obj_ordered@colData)

write.csv(partitioned_obj_ordered@colData$pseudotime, file.path(figs_dir_path, paste0("metadata_pseudotime_clus.csv")))


## Save objects
# save_monocle_objects(cds=partitioned_obj_ordered, 
#                     directory_path=file.path(figs_dir_path, paste0(filename,"_root_mixedcycling")), 
#                     comment='GBM_GTE_pseudotime root by cluster, 20250320')
saveRDS(partitioned_obj_ordered, file=file.path(parent_dir_path_obj, paste0(filename,"_root_mixedcycling.rds")))

#### ===================================================================== ####
#### Compute Pseudotime: Select root cells by function ####
#### ===================================================================== #### 
## FUNCTION: programmatically choose which cells come first in the trajectory
get_earliest_principal_node <- function(cds_obj, column, cluster){
  cell_ids <- which(colData(cds_obj)[, column] == cluster)
  closest_vertex <- cds_obj@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds_obj), ])
  root_pr_nodes <- igraph::V(principal_graph(cds_obj)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}

## Order Cells
partitioned_obj_autoroot <- order_cells(partitioned_obj, root_pr_nodes=get_earliest_principal_node(partitioned_obj, cluster_col, "Mixed Cycling"))

## Export Pseudotime graph
p <- partitioned_obj_autoroot %>% plot_cells(
  color_cells_by = "pseudotime", 
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE,
  label_roots = FALSE,
  trajectory_graph_color = "black"
  # graph_label_size = 5,
)
ggsave(file.path(figs_dir_path, paste0(filename, "_pseuodtime_autoroot.tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

p <- partitioned_obj_autoroot %>% plot_cells(
  color_cells_by = "pseudotime", 
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE,
  label_roots = TRUE,
  trajectory_graph_color = "black",
  graph_label_size = 3
)
ggsave(file.path(figs_dir_path, paste0(filename, "_pseuodtime_autoroot_labelroot.tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

p <- partitioned_obj_autoroot %>% plot_cells(
  color_cells_by = "pseudotime", 
  label_cell_groups = FALSE, 
  label_leaves = TRUE, 
  label_branch_points = FALSE,
  label_roots = FALSE,
  trajectory_graph_color = "black",
  graph_label_size = 3
)
ggsave(file.path(figs_dir_path, paste0(filename, "_pseuodtime_autoroot_labelleaves.tiff")), 
  plot = p, units="in", width=size*1.2, height=size, dpi=300, compression = 'lzw')

# save pseudotime ordering
pseudotime <- pseudotime(partitioned_obj_autoroot) 	
partitioned_obj_autoroot@colData$pseudotime_auto <- pseudotime

glimpse(partitioned_obj_autoroot@colData)

write.csv(partitioned_obj_autoroot@colData$pseudotime_auto, file.path(figs_dir_path, paste0("metadata_pseudotime_auto.csv")))

## Save objects
# save_monocle_objects(cds=partitioned_obj_autoroot, 
#                     directory_path=file.path(figs_dir_path, paste0(filename,"_autoroot")), 
#                     comment='GBM_GTE_pseudotime auto choose root cell, 20250320')
saveRDS(partitioned_obj_autoroot, file=file.path(parent_dir_path_obj, paste0(filename,"_autoroot.rds")))

### End of Script
sessionInfo()
