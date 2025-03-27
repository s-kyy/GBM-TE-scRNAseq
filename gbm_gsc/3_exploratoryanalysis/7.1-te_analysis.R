#!usr/bin/env Rscript

library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("At least 2 filepaths must be supplied: [gte.rds] [te_genes] [res] [figure_path]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2]) ) { 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_name_obj <- basename(parent_dir_path_obj)
    te_path <- args[2] # e.g. "GRCh38_Ensembl_rmsk_TE_v23.gtf"
  } else {
    stop("Filepaths provided do not exist. Exiting...", call.=FALSE)
  }
  
  if (length(args)>=3) {
    # Name of column from integrated dataset 
    cluster_res_num <- as.numeric(args[3])
    cluster_res <- paste0("integrated_snn_res.", cluster_res_num)
    # Custom column name for celltypes in metadata
    cluster_res_symbol <- unlist(strsplit(as.character(cluster_res_num), ".", fixed=TRUE))
    cluster_col <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_celltypes")
    if (grepl("gbm", parent_dir_path_obj, fixed = TRUE)) {
      cluster_col_gsc <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_gsctypes")
      cluster_col_gsc_cc <- paste0(cluster_col_gsc,"_cc")
      cluster_col_gbm_neftel <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_gbmneftel")
      cluster_col_gbm_tcga <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_gbmtcga")
    }
  } else {
    cluster_res <- "integrated_snn_res.0.4"
    cluster_col <- "int04_celltypes"
    if (grepl("gbm", parent_dir_path_obj, fixed = TRUE)) {
      cluster_col_gsc <- "int04_gsctypes"
      cluster_col_gsc_cc <- "int04_gsctypes_cc"
      cluster_col_gbm_neftel <- "int04_gbmneftel"
      cluster_col_gbm_tcga <- "int04_gbmtcga"
    }
  }

  if (length(args)>=4) {
    figs_dir_name <- args[4]
  } else {
    figs_dir_name <- "figs_te_analysis"
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
library(viridis)

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(108)

#### ===================================================================== ####
#### Load Datasets ####
#### ===================================================================== ####
seurat_obj <- readRDS(obj_path) 
DefaultAssay(seurat_obj) <- "RNA"
# cluster_markers <- read.csv(marker_path)
size    <- 5

# Import TE genes
te_df <- read.table(te_path, sep=",", header = TRUE,
  col.names=c("gene_id", "class", "family"), fill=FALSE, strip.white=TRUE)
te_df$gene_id <- sub("_", "-", te_df$gene_id) # in seurat object, `_` were replaced with `-`
rownames(te_df) <- te_df$gene_id
dim(te_df)
head(te_df)

# Categorize broad Retrotransposon types
te_df$type <- "Other"
te_df$type[which((te_df$class == "LINE") )] <- "LINE"
te_df$type[which((te_df$class == "SINE") | (te_df$class == "Retroposon") )] <- "SINE" #Alu and SVA elements
te_df$type[which((te_df$class == "LTR") )] <- "LTR"

table(te_df$class)
table(te_df$type)
# LINE   LTR Other  SINE 
#  171   578   374    57

print("SVA elements")
te_df$gene_id[which(te_df$class == "Retroposon")]
print("SINE elements")
te_df$gene_id[which(te_df$class == "SINE")]
print("LINE elements")
te_df$gene_id[which(te_df$type == "LINE")]
print("LTR elements")
te_df$gene_id[which(te_df$type == "LTR")]
print("Other elements")
te_df$gene_id[which(te_df$type == "Other")]

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

figs_dir_path <- file.path(parent_dir_path_obj, figs_dir_name)

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

#### ===================================================================== ####
#### % of Retrotransposons ####
#### ===================================================================== ####

# Number of RTs detected in dataset
matched <- intersect(te_df$gene_id, rownames(seurat_obj$RNA@data))
non.matched <- te_df$gene_id[!te_df$gene_id %in% matched]
print(paste0("Number of detected TEs: ",length(matched)))
print(paste0("Number of undetected TEs: ",length(non.matched)))

table(te_df$type[which(te_df$gene_id %in% matched)])
# LINE   LTR Other  SINE
#  161   565   308    57

te_df <- filter(te_df, gene_id %in% matched)
total_num_te <- length(te_df$gene_id)

meta <- seurat_obj@meta.data

# Count UMI of TE-types per gsctypes (counts)
dim(seurat_obj)
dim(meta)

summarize_te_distributions <- function(seurat_obj_, cluster_var, te_df_) {
  summary_counts <- unique(seurat_obj_[[cluster_var]])
  rownames(summary_counts) <- summary_counts[,1]
  summary_counts$LINE <- NA
  summary_counts$SINE <- NA
  summary_counts$LTR <- NA
  summary_counts$Other <- NA 

  summary_data <- unique(seurat_obj_[[cluster_var]])
  rownames(summary_data) <- summary_data[,1]
  summary_data$LINE <- NA
  summary_data$SINE <- NA
  summary_data$LTR <- NA
  summary_data$Other <- NA 

  summary_numte <- unique(seurat_obj_[[cluster_var]])
  rownames(summary_numte) <- summary_numte[,1]
  summary_numte$LINE <- NA
  summary_numte$SINE <- NA
  summary_numte$LTR <- NA
  summary_numte$Other <- NA 

  summary_te_celltypes <- unique(seurat_obj_[[cluster_var]])
  rownames(summary_numte) <- summary_numte[,1]
  summary_numte$LINE <- NA
  summary_numte$SINE <- NA
  summary_numte$LTR <- NA
  summary_numte$Other <- NA 

  summary_cells <- te_df_
  summary_cells$num_celltypes <- 0
  # summary_cells$num_cells_expressing_te <- 0

  for ( cluster in levels(meta[[cluster_var]]) ) {
    
    barcodes <- rownames(filter(meta, !!sym(cluster_var) == cluster))

    for ( te_type in c("LINE", "SINE", "LTR", "Other") ) {
      
      print(paste("Cluster: ", cluster, ", TE type: ", te_type))
      te_subset <- te_df_$gene_id[which(te_df_$type==te_type)]
      te_subset_size <- length(te_df_$gene_id[which(te_df_$type==te_type)])

      # sum of raw counts across TE subset and cluster
      summary_counts[cluster,te_type] <- sum(Matrix::rowSums(seurat_obj_[te_subset,barcodes]$RNA@counts))
      
      # sum of log-normalized counts across TE subset and cluster
      summary_data[cluster,te_type] <- sum(Matrix::rowSums(seurat_obj_[te_subset,barcodes]$RNA@data))
      
      # total number of TE in subset - num TE not expressed in cluster/celltype
      summary_numte[cluster,te_type] <- te_subset_size - sum(tabulate(seurat_obj_[te_subset,barcodes]$RNA@counts@i + 1) == 0) 
    }
    
    print(paste("Number of clusters expression TE"))
    summary_cells[[cluster]] <- 0
    
    for ( te_id in summary_cells$gene_id ) {
      
      num_cells_expressing_te <- sum(Matrix::colSums(seurat_obj_[te_id,barcodes]$RNA@counts) >= 1)
      
      if ( num_cells_expressing_te > 0 ) { 
        summary_cells[te_id,"num_celltypes"] <- summary_cells[te_id,"num_celltypes"] + 1
        summary_cells[te_id,cluster] <- num_cells_expressing_te
      }
    }
  }

  write.csv(summary_counts, file=file.path(figs_dir_path,paste0("summary_te_counts_",cluster_var,".csv")), row.names = FALSE)
  write.csv(summary_data, file=file.path(figs_dir_path,paste0("summary_te_data_",cluster_var,".csv")), row.names = FALSE)
  write.csv(summary_numte, file=file.path(figs_dir_path,paste0("summary_te_numte_",cluster_var,".csv")), row.names = FALSE)
  write.csv(summary_cells, file=file.path(figs_dir_path,paste0("te_df_numclusters_",cluster_var,".csv")), row.names = FALSE)
  return (TRUE)
}

if (grepl("gbm", parent_dir_path_obj, fixed = TRUE)) {
  # summarize_te_distributions(seurat_obj, cluster_col_gsc, te_df)
  # summarize_te_distributions(seurat_obj, cluster_col, te_df)
  summarize_te_distributions(seurat_obj, cluster_col_gbm_neftel, te_df)
  summarize_te_distributions(seurat_obj, cluster_col_gbm_tcga, te_df)
} else {
  summarize_te_distributions(seurat_obj, cluster_col, te_df)
}


### End of Script
sessionInfo()
