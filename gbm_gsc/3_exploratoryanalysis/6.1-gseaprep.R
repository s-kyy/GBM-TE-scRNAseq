#!usr/bin/env Rscript

.libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 2 filepath must be supplied: [markers.csv] [fig_dir_name]", call.=FALSE)
} else if (length(args)>=1) {
  
  # verify filepaths
  if (file.exists(args[1])) { 
    marker_path <- args[1]
    filename <- basename(path_ext_remove(marker_path))
  } else {
    stop("Filepaths provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>1) {
    figs_dir_name <- args[2]
  } else {
    figs_dir_name <- "gseaprep"
  }

} else {
    stop("Error in calling R or filepaths. Exiting...", call.=FALSE)
}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
set.seed(108)

# library(Matrix)
library(tidyverse) 
library(dplyr)
# library(rlang)
library(stringr)

set.seed(108)
options(warn=1) #print warning messages as they occur

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
markers_df <- read.csv(marker_path)
colnames_list <- unique(markers_df$cluster_col)

subdir <- file.path(getwd(), paste0(format(Sys.Date(), "%Y%m%d"), "_", filename,"_avgexp_bygroup"))

ifelse(!dir.exists(file.path(subdir)),
        dir.create(file.path(subdir),recursive=T),
        "Directory Exists")

figs_dir_path <- file.path(subdir,figs_dir_name)

ifelse(!dir.exists(file.path(figs_dir_path)),
        dir.create(file.path(figs_dir_path),recursive=T),
        "Directory Exists")

#### =========================================== ####
#### A. Prepare Input for GSEA: expression tables
#### =========================================== ####

# extract gene names from dataframe into new column
if ("X" %in% colnames(markers_df)) {
  markers_df$gene <- sapply(strsplit(
                            markers_df$X,split = "...", fixed = TRUE), function(x) x[1] )
  markers_df$X <- NULL
} else if (rownames(markers_df)[1] != "1") {
  markers_df$gene <- rownames(markers_df)
}

for (idents in colnames_list) {
    
  print("Ranking genes ...")

  cluster_list <-  unique( as.character( markers_df$cluster[which(markers_df$cluster_col == idents)] ) )
  
  for ( group in cluster_list ) {
    # Create list of ranked genes 
    # Rank genes : log2FC * -log10(p+ 1e-310) (with unadjusted pvalues)
    ranked_genes <- markers_df %>% 
      dplyr::select(gene,avg_log2FC,p_val,cluster,cluster_col) %>%        # filter columns log2FC, p-val, genes, cluster, colname
      dplyr::filter(cluster_col == idents & cluster == group &            # filter genes (log2FC >0.1, p-val (non-adjusted) < 0.01)
                    p_val <= 0.01 & (avg_log2FC >= 0.1 | avg_log2FC <= -0.1)) %>% 
      dplyr::mutate( rank = avg_log2FC * -(log10(p_val + 1e-310)) ) %>%   # Calculate ranking metric
      dplyr::arrange(desc(rank))                                          # Sort genes by ranking metric
    print(head(ranked_genes))
    print("Exporting Ranked genes ...")
    group_name <- str_replace_all(group,"/","-")
    write.csv(ranked_genes, file.path(figs_dir_path, paste0(filename,"_rankedgenes_",idents,"_",group_name,".csv")), row.names = FALSE)

    ranked_genes_input <- ranked_genes %>% dplyr::select(gene,rank)
    write.table(ranked_genes_input, file.path(figs_dir_path, paste0(filename,"_rankedgenes_",idents,"_",group_name,".rnk")),
       sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

    print("Ranked genes exported...")
    
  }
}
    
## Session Info
sessionInfo()
