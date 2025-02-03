#!usr/bin/env Rscript

# Inputs
# Counts matrix 
# cell annotation file = tab-delimited, no headers.
# Columns: cell barcode, celltype
# gene annotation file = tab-delimited, no headers. See `4.1-ordergenes.sh`
# Columns: gene name, chromosome, and gene span
# https://github.com/broadinstitute/inferCNV/wiki/File-Definitions

#### Parse Arguments ####
.libPaths(c("~/scratch/tcga-gbm-R4.1-lib/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()

library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("At least 4 args must be supplied: [xxx.rds] [gene_annotations.txt] [celltype_col] [num_threads/workers] ", call.=FALSE)
} else {
  
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2])){ 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_name_obj <- basename(parent_dir_path_obj)
    gene_annotations_path <- args[2]
  } else {
    stop("Filepath provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>2) {
    celltype_col <- args[3]
  } else {
    celltype_col <- "int06_celltype"
  }

  if (length(args)>3) {
    threads <- as.numeric(args[4])
  } else {
    threads <- 6
  }
  

}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
set.seed(108)

library(Matrix)
library(ggplot2) # BiocManager (Bioconductor version 3.14)
library(tidyverse) # BiocManager (Bioconductor version 3.14)
library(infercnv) # R 4.1
library(future)

options(warn=1) # print warning messages as they occur

set.seed(108, kind = "L'Ecuyer-CMRG") 

### To avoid this error:
# Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) : 
# The total size of the x globals exported for future expression ('FUN()') is yyy GiB.. This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
## Set max size to 10.0 GB
options(future.globals.maxSize= 10000 * 1024^2)  # 1000 * 1024^2 = 1Gb
options(future.rng.onMisuse = "ignore")

future::plan(sequential)
print(plan())


#### =========================================== ####
#### Load Objects ####
#### =========================================== ####

seurat_obj <- readRDS(obj_path) 

subdir <- file.path(parent_dir_path_obj,"infercnv")

ifelse(!dir.exists(file.path(subdir)),
        dir.create(file.path(subdir),recursive=T),
        "Directory Exists")

# generate and export cell annotations file
cell_annotations_path <- file.path(subdir,"cell_annotations.txt")
cell_annotations <- seurat_obj@meta.data[celltype_col]
write.table(cell_annotations, cell_annotations_path,
  sep="\t", col.names=FALSE, quote=FALSE, row.names=TRUE
)

#### =========================================== ####
#### Infer CNVs ####
#### =========================================== ####

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix= seurat_obj@assays$RNA@counts[,colnames(seurat_obj)],
                      # GetAssayData(seurat_obj, slot="counts", assay="RNA"),
  annotations_file= cell_annotations_path,
  gene_order_file=gene_annotations_path,
  delim="\t",
  ref_group_names=c("Immune","Oligodendrocyte", "Endothelia") #'Microglia' (removed as reference after preliminary run)
) 

future::plan(multisession, workers = threads)
print(plan())

infercnv_obj <- infercnv::run(
  infercnv_obj,
  out_dir=subdir, 
  cutoff=0.1,  # cutoff=1 (Smart-seq2), cutoff=0.1 (10x Genomics)
  window_length=101, #default=101
  min_cells_per_gene=10, #Bhaduri et al. 2020 
  cluster_by_groups=TRUE, 
  denoise=TRUE,
  HMM=TRUE, HMM_type="i6", 
  resume_mode=FALSE,
  png_res=300, #dpi
  num_threads=threads #default=4
)

saveRDS(infercnv_obj, file=file.path(subdir, paste0(filename,"_infercnv.rds")))

#### End of Script
sessionInfo()