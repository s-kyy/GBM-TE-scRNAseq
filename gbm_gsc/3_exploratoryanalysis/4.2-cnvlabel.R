#!usr/bin/env Rscript

# .libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
# .libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least 2 filepaths must be supplied: [xxx_ge.rds] [yyy_gte.rds] [zzz-sucblusters.cell_groupings] [figs_dir_name] [cluster_res]", call.=FALSE)
} else if (length(args)>=2) {
  
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2]) && file.exists(args[3])) { 
    obj_path <- args[1] 
    obj_path_2 <- args[2] 
    groupings_path <- args[3] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_path_obj_2 <- dirname(obj_path_2)
    # parent_dir_name_obj <- basename(parent_dir_path_obj)
  } else {
    stop("Filepaths provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>3) {
    figs_dir_name <- args[4]
  } else {
    figs_dir_name <- "figs_celltypes_annotated"
  }

  if (length(args)==5) {
    # Name of column from integrated dataset 
    cluster_res_num <- as.numeric(args[5])
    cluster_res <- paste0("integrated_snn_res.", cluster_res_num)
    # Custom column name for celltypes in metadata
    cluster_res_symbol <- unlist(strsplit(as.character(cluster_res_num), ".", fixed=TRUE))
    cluster_col <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_celltypes")
    cluster_col_gsc <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_gsctypes")
    cluster_col_cnv <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_tumour")
  } else {
    cluster_res <- "integrated_snn_res.0.4"
    cluster_col <- "int04_celltypes"
    cluster_col_gsc <- "int04_gsctypes"
    cluster_col_cnv <- paste0("int04_tumour")
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
library(dplyr)
library(viridis)
library(rlang)
library(reshape2)
library(car)

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(108)
options(warn=1) #print warning messages as they occur

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
seurat_obj_ge <- readRDS(obj_path) 
seurat_obj_gte <- readRDS(obj_path_2) 
cellgroupings <- read.table(groupings_path, sep="\t", header=TRUE)
row.names(cellgroupings) <- cellgroupings[,2]
cellgroupings <- cellgroupings %>% separate(col= 1, into=c("maingroup", "subcluster"), sep="\\.")
head(cellgroupings)
size    <- 5
ndims <- 25

# subdir <- file.path(getwd(), paste0(format(Sys.Date(), "%Y%m%d"), "_", sample_name,"_", filename))
subdir <- parent_dir_path_obj
subdir_2 <- parent_dir_path_obj_2

ifelse(!dir.exists(file.path(subdir)),
        dir.create(file.path(subdir),recursive=T),
        "Directory Exists")
ifelse(!dir.exists(file.path(subdir_2)),
        dir.create(file.path(subdir_2),recursive=T),
        "Directory Exists")

figs_dir_path <- file.path(subdir, paste0(figs_dir_name, "_ge"))
figs_dir_path_2 <- file.path(subdir_2, paste0(figs_dir_name, "_gte"))

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")
ifelse(!dir.exists(figs_dir_path_2),
        dir.create(figs_dir_path_2,recursive=T),
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
#### Label Tumour cells ####
#### =========================================== ####

# Tumour cell labels
nontumour_cells <- c(
  "SF11285_s2",
  "Immune",
  "Oligodendrocyte", 
  "Endothelia", 
  "Microglia"
)

# create new column to label tumour cells
meta <- seurat_obj_ge@meta.data
meta[[cluster_col_cnv]] <- 1

glimpse(meta)

# label non-tumour cells belonging to one of the nontumour_cells groups
for (cell in cellgroupings[,"cell"]){
  if ((cellgroupings[cell,"maingroup"] %in% nontumour_cells) || (cellgroupings[cell,"subcluster"] %in% nontumour_cells)) {
    meta[cell,cluster_col_cnv] <- 0
  } 
}

seurat_obj_ge <- AddMetaData(seurat_obj_ge, meta[cluster_col_cnv], col.name = cluster_col_cnv)
seurat_obj_gte <- AddMetaData(seurat_obj_gte, meta[cluster_col_cnv], col.name = cluster_col_cnv)

glimpse(seurat_obj_ge@meta.data)
glimpse(seurat_obj_gte@meta.data)

if (class(meta$mitoRatio) == "data.frame") { 
  # meta$mitoRatio <- meta$mitoRatio[,1] 
  meta$mitoRatio <- NULL
} else if (class(meta$mitoRatio) == "list") {
  # meta$mitoRatio <- meta$mitoRatio[[1]]
  meta$mitoRatio <- NULL
}
meta <- meta %>% dplyr::mutate(across(everything(), as.character))
write.csv(meta, file.path(figs_dir_path, paste0("metadata_ge_cnv.csv")))

meta <- seurat_obj_gte@meta.data
if (class(meta$mitoRatio) == "data.frame") { 
  # meta$mitoRatio <- meta$mitoRatio[,1] 
  meta$mitoRatio <- NULL
} else if (class(meta$mitoRatio) == "list") {
  # meta$mitoRatio <- meta$mitoRatio[[1]]
  meta$mitoRatio <- NULL
}
meta <- meta %>% dplyr::mutate(across(everything(), as.character))
write.csv(meta, file.path(figs_dir_path_2, paste0("metadata_gte_cnv.csv")))
rm("meta")

if (class(seurat_obj_ge@meta.data$mitoRatio) == "data.frame" || class(seurat_obj_ge@meta.data$mitoRatio) == "list") { 
  seurat_obj_ge@meta.data$mitoRatio <- NULL 
} 
if (class(seurat_obj_gte@meta.data$mitoRatio) == "data.frame"|| class(seurat_obj_gte@meta.data$mitoRatio) == "list") { 
  seurat_obj_gte@meta.data$mitoRatio <- NULL 
} 

saveRDS(seurat_obj_ge, file = file.path(subdir, paste0("gbm_ge_celltypes_cnv.rds")))
saveRDS(seurat_obj_gte, file = file.path(subdir_2, paste0("gbm_gte_celltypes_cnv.rds")))

#### ===================================================================== ####
#### UMAP plots by tumour label ####
#### ===================================================================== ####

barcode_cnv <- rownames(filter(seurat_obj_ge@meta.data, !!sym(cluster_col_cnv) == 1))

p <- DimPlot(seurat_obj_ge, reduction = "umap", 
            cells.highlight = WhichCells(object = seurat_obj_ge, cells = barcode_cnv) )
ggsave(file.path(figs_dir_path, paste0("tumourcells",cluster_col_cnv,"_UMAP.tiff")),
      plot = p, units="in", width=size*1.2, height=size*1.2, dpi=300, compression = 'lzw')
p <- DimPlot(seurat_obj_gte, reduction = "umap", 
            cells.highlight = WhichCells(object = seurat_obj_gte, cells = barcode_cnv) )
ggsave(file.path(figs_dir_path_2, paste0("tumourcells",cluster_col_cnv,"_UMAP.tiff")),
      plot = p, units="in", width=size*1.2, height=size*1.2, dpi=300, compression = 'lzw')

print("UMAPs by tumour cell label exported")

#### ===================================================================== ####
#### CSV files ####
#### ===================================================================== ####

calcProportions <- function (meta, col1, col2) {
  grouped_counts_table <- table(meta[c(col1, col2)])
  grouped_counts_table <- cbind(grouped_counts_table, prop.table(grouped_counts_table, margin=1)*100)
  return(grouped_counts_table)
}

write.csv(calcProportions(seurat_obj_ge@meta.data, "sample",cluster_col_cnv), file.path(figs_dir_path, paste0("summary_sampleVS",cluster_col_cnv,".csv")), row.names=FALSE)
write.csv(calcProportions(seurat_obj_ge@meta.data, cluster_col,cluster_col_cnv), file.path(figs_dir_path, paste0("summary_",cluster_col,"VS",cluster_col_cnv,".csv")), row.names=FALSE)
write.csv(calcProportions(seurat_obj_ge@meta.data, cluster_col_gsc,cluster_col_cnv), file.path(figs_dir_path, paste0("summary_",cluster_col_gsc,"VS",cluster_col_cnv,".csv")), row.names=FALSE)
# write.csv(table(seurat_obj_ge@meta.data["sample",cluster_col,cluster_col_cnv]), file.path(figs_dir_path, paste0("summary_sampleVS",cluster_col_cnv,"VS",cluster_col,".csv")), row.names=FALSE)
# write.csv(table(seurat_obj_ge@meta.data["sample",cluster_col_gsc,cluster_col_cnv]), file.path(figs_dir_path, paste0("summary_sampleVS",cluster_col_cnv,"VS",cluster_col_gsc,".csv")), row.names=FALSE)
# write.csv(table(seurat_obj_ge@meta.data["sample",cluster_col_gsc_cc,cluster_col_cnv]), file.path(figs_dir_path, paste0("summary_sampleVS",cluster_col_cnv,"VS",cluster_col_gsc_cc,".csv")), row.names=FALSE)

### End of Script
sessionInfo()
