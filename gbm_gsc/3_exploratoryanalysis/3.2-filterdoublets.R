#!usr/bin/env Rscript

#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 1 filepath must be supplied: [xxx.rds] [sample_name]", call.=FALSE)
} else if (length(args)<2) {
  # verify filepaths
  if (file.exists(args[1])){ 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    sample_name <- ""    
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_name_obj <- basename(parent_dir_path_obj)
  } else {
    stop("Filepath provided does not exist. Exiting...", call.=FALSE)
  }
} else if (length(args)<=2) {
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    sample_name <- args[2]
}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
.libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()

set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 
library(DoubletFinder)
# DoubletFinder@03e9f37f891ef76a23cc55ea69f940c536ae8f9f (April 10, 2024)

set.seed(108)

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
# Ensure objects have reached the RunUMAP step
seurat.obj <- readRDS(obj_path) 
size    <- 5

figs_dir_path <- file.path(parent_dir_path_obj, "figs_celltypes")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

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
print(paste("Minimum PC with a difference in variation of 10% compared to next PC:", pc_10))
co2

# Minimum of the two calculation
min_pc <- min(pc_most_var, pc_10)
print(paste("Minumum PC between the two options:", min_pc))

plot_df <- data.frame(dimensions = 1:length(pct_var_per_pc),
           stdev = seurat.obj[["pca"]]@stdev,
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

#### =========================================== ####
#### Split objects ####
#### =========================================== ####

if (grepl("healthy", parent_dir_name_obj, fixed = TRUE)) {
  # split healthy dataset by sample. 
  obj_list <- SplitObject(object = seurat.obj, split.by="sample")
} else if (grepl("gbm", parent_dir_name_obj, fixed = TRUE)) {
  # copy GBM_WANG2020 sample names to sample_orig column
  for (i in which(is.na(seurat.obj@meta.data$sample_orig))){
    seurat.obj@meta.data$sample_orig[i] <- as.character(seurat.obj@meta.data$sample[i])
  }
  # split gbm_bhaduri sample by lane and gbm_wang by sample
  obj_list <- SplitObject(object = seurat.obj, split.by="sample_orig")
}

#### =========================================== ####
#### Run DoubletFinder ####
#### =========================================== ####

findDoublets <- function(obj, dims) {
  #### pK Identification (no ground-truth) 
  sweep.res.list <- paramSweep(obj, PCs = 1:min_pc, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats)

  #### Homotypic Doublet Proportion Estimate 
  homotypic.prop <- modelHomotypic(annotations)           
    ## ex: annotations <- obj@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(obj@meta.data))  
    ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  #### Run DoubletFinder with varying classification stringencies 
  result <- doubletFinder(obj, PCs = 1:min_pc, 
                          pN = 0.25, pK = 0.09, nExp = nExp_poi, 
                          reuse.pANN = FALSE, sct = FALSE)
  result <- doubletFinder(obj, PCs = 1:min_pc, 
                          pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, 
                          reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
  return(result)
}

for (i in 1:length(obj_list)) {
  obj_list[[i]] <- findDoublets(obj_list[[i]], min_pc)
}
prin("Doublets annotated")

# merge dataset
prin("Seurat Object merged")

# export object
prin("Seurat Object exported")


# create UMAP with doublets labeled
prin("Doublets annotated in UMAP")

# export csv with cell counts per sample and number of labeled doublets
print("csv exported")


#### End of Script ####
sessionInfo()