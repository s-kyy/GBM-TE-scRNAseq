#!usr/bin/env Rscript

.libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 2 filepaths must be supplied: [xxx.rds] [yyy.rds]", call.=FALSE)
} else if (length(args)>=2) {
  
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2])) { 
    obj_path <- args[1] 
    obj_path_2 <- args[2] 
    filename <- basename(path_ext_remove(obj_path))
    filename_2 <- basename(path_ext_remove(obj_path_2))
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_path_obj_2 <- dirname(obj_path_2)
    # parent_dir_name_obj <- basename(parent_dir_path_obj)
  } else {
    stop("Filepaths provided does not exist. Exiting...", call.=FALSE)
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
seurat_obj_ge <- readRDS(obj_path) 
seurat_obj_gte <- readRDS(obj_path_2)
size    <- 5
ndims <- 25

# subdir <- file.path(getwd(), paste0(format(Sys.Date(), "%Y%m%d"), "_", sample_name,"_", filename))

figs_dir_path <- file.path(parent_dir_path_obj, "figs_postfilterDeadCells")
figs_dir_path_2 <- file.path(parent_dir_path_obj_2, "figs_postfilterDeadCells")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")
ifelse(!dir.exists(figs_dir_path_2),
        dir.create(figs_dir_path_2,recursive=T),
        "Directory Exists")

sample_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", 
  "#F0E442", "#CC79A7", "#ff716e",
  "#999999", "#0072B2", "#194c76")
female_samples <- "SF11209"

## ========================================= ##
## Filter out Doublets ##
## ========================================= ##

DefaultAssay(seurat_obj_ge) <- "RNA"
DefaultAssay(seurat_obj_gte) <- "RNA"

ncells_old <- dim(seurat_obj_ge)[2]
seurat_obj_ge <- subset(
  x = seurat_obj_ge, 
  subset = (integrated_snn_res.0.4 != 12)
) 
print(paste(ncells_old - dim(seurat_obj_ge)[2], "cells were removed from gbm ge object"))

barcode <- colnames(seurat_obj_ge)
seurat_obj_gte <- subset(seurat_obj_gte, cells=barcode)
print(paste(ncells_old - dim(seurat_obj_gte)[2], "cells were removed from gbm gte object"))

saveRDS(seurat_obj_ge, file = file.path(parent_dir_path_obj, paste0(filename, "_filtDC.rds")))
saveRDS(seurat_obj_gte, file = file.path(parent_dir_path_obj_2, paste0(filename_2, "_filtDC.rds")))

#### =========================================== ####
#### Make BarPlot (post filter-DeadCells) ####
#### =========================================== ####
makePlots <- function(obj, col, file_name, output_path) {
  print(paste("Making figures and summary csvs for:", file_name))
  DefaultAssay(obj) <- "RNA"

  # Factor sample names
  if (is.factor(obj@meta.data$orig.ident)) {
    obj@meta.data$orig.ident <- droplevels(obj@meta.data$orig.ident)
  }
  if (is.factor(obj@meta.data$sample)) {
    sample_names <- sort(unique(droplevels(obj@meta.data$sample)))
  } else {
    sample_names <- sort(unique(obj@meta.data$sample))
  }
  obj@meta.data$sample <- factor(
    obj@meta.data$sample, 
    levels = sample_names
  )
  sample_names_indexed <- sample_names
  names(sample_names_indexed) <- seq(1,length(sample_names))

  print("Bar plot Cell Count per Sample")
  p <- obj@meta.data %>%
      ggplot(aes(x=sample, fill=sample)) + 
      labs(x = "", y = "Cell Count") +
      scale_fill_manual(values = col) + 
      labs(fill="Sample ID") + 
      geom_bar() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=10), legend.position = "none")
  ggsave(file.path(output_path, paste0("cellcounts_filtDeadCells.tiff")), 
        plot = p, units="in", width=size*1, height=size*0.8, dpi=300, compression = 'lzw')

  # Table of Cell Counts per Sample
  df <- obj@meta.data %>% 
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
  obj.list <- SplitObject(obj, split.by="sample")
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
  write.csv(df, file= file.path(output_path, paste0("samplecounts_filtDeadCells.csv")), row.names = F)
}

makePlots(obj = seurat_obj_ge, col=sample_palette, file_name=filename, output_path=figs_dir_path)
makePlots(obj = seurat_obj_gte, col=sample_palette, file_name=filename_2, output_path=figs_dir_path_2)

#### End of Script ####
sessionInfo()
