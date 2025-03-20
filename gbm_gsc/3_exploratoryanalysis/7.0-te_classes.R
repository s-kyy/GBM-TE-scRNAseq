#!usr/bin/env Rscript

library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 3 filepaths must be supplied: [te_genes.gtf] [figure_path]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1])) { 
    gtf_path <- args[1] # e.g. "GRCh38_Ensembl_rmsk_TE_v23.gtf"
  } else {
    stop("Filepaths provided do not exist. Exiting...", call.=FALSE)
  }
  
  if (length(args)>=2) {
    figs_dir_name <- args[2]
  } else {
    figs_dir_name <- "meta"
  }
} 

# Import Packages
library(rtracklayer)
library(data.table)

# Set output path
figs_dir_path <- file.path(getwd(), figs_dir_name)

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

## Obtain list of Retrotransposon genes
gtf_gr <- rtracklayer::import(gtf_path) # creates a GRanges object
gtf_df <- as.data.frame(gtf_gr)
genes <- unique(gtf_df[ ,c("gene_id", "class_id", "family_id")])
## save retrotransposon genes in a text file
fwrite(genes, file=file.path(figs_dir_path, paste0(basename(path_ext_remove(gtf_path)),"_geneNames.txt")), sep=",")

### End of Script
sessionInfo()
