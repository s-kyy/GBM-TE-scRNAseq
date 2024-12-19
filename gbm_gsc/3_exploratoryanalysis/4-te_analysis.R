
# import 2 files, seurat object + cellular marker CSV (one at a time)

library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("At least 3 filepaths must be supplied: [xxx.rds] [xxx_markers_0.3.csv] [xxx_markers_0.4.csv] [te_genes.gtf]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1]) & file.exists(args[2]) & file.exists(args[3])) { 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_name_obj <- basename(parent_dir_path_obj)
    marker_path3 <- args[2]
    marker_path4 <-args[3]
#     parent_dir_path_marker <- dirname(marker_path)
#     parent_dir_name_marker <- basename(parent_dir_path_marker)
  } else {
    stop("Filepaths provided do not exist. Exiting...", call.=FALSE)
  }
}

#### ===================================================================== ####
#### Import Packages ####
#### ===================================================================== ####

set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 
library(dplyr)
library(reshape2)
library(viridis)
library(ggrepel)

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(108)

#### ===================================================================== ####
#### Load Datasets ####
#### ===================================================================== ####
seurat.obj <- readRDS(obj_path) 
cluster_markers3 <- read.csv(marker_path3)
cluster_markers4 <- read.csv(marker_path4)
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
  avg_exp_threshold <- 10

} else if (grepl("gbm", parent_dir_name_obj, fixed = TRUE)) {
  sample_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#CC79A7", "#ff716e",
    "#999999", "#0072B2", "#194c76")
  female_samples <- "SF11209"
  avg_exp_threshold <- 25
}

seurat.obj$`integrated_snn_res.0.3` <- factor(seurat.obj$`integrated_snn_res.0.3`, levels= sort(unique(seurat.obj$`integrated_snn_res.0.3`)))
seurat.obj$`integrated_snn_res.0.4` <- factor(seurat.obj$`integrated_snn_res.0.4`, levels= sort(unique(seurat.obj$`integrated_snn_res.0.4`)))

figs_dir_path <- file.path(parent_dir_path_obj, "figs")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

print(seurat.obj)
