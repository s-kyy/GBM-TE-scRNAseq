#!usr/bin/env Rscript

library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least 3 filepaths must be supplied: [cluster_markers.csv] [te_df] [figure_path] [celltypes]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2]) ) { 
    marker_path <- args[1] 
    filename <- basename(path_ext_remove(marker_path))
    parent_dir_path_obj <- dirname(marker_path)
    parent_dir_name_obj <- basename(parent_dir_path_obj)
    # marker_path <- args[2]
    te_path <- args[2] # e.g. "GRCh38_Ensembl_rmsk_TE_v23.gtf"
  } else {
    stop("Filepaths provided do not exist. Exiting...", call.=FALSE)
  }
  
  if (length(args)>=3) {
    figs_dir_name <- args[3]
  } else {
    figs_dir_name <- "figs_te_analysis"
  }
  if (length(args)>=4) {
    cluster_group <- args[4]
  } else {
    cluster_group <- "celltypes"
  }
}


#### ===================================================================== ####
#### Import Packages ####
#### ===================================================================== ####

set.seed(108)

library(ggplot2)
library(tidyverse) 
# library(reshape2)
library(dplyr)
library(viridis)
library(ggrepel)

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(108)

#### ===================================================================== ####
#### Load Datasets ####
#### ===================================================================== ####
cluster_markers <- read.csv(marker_path)
te_df <- read.csv(te_path)
te_df <- te_df[,1:4]
te_df$type <- factor(te_df$type, levels=c("LINE", "SINE", "LTR", "Other"))
size    <- 5

if (grepl("healthy", parent_dir_path_obj, fixed = TRUE)) {
  sample_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#CC79A7", "#ff716e",
    "#999999", "#0072B2", "#194c76", 
    "#D55E00", "#3a4f41", "#6699cc", "#713e5a")
  female_samples <- c("SRR9262922", "SRR9262937",
                        "SRR9264382", "SRR9264383",
                        "SRR9264388")
  avg_exp_threshold <- 10

} else if (grepl("gbm", parent_dir_path_obj, fixed = TRUE)) {
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
#### Relabel Celltypes ####
#### ===================================================================== ####

if (grepl("healthy", parent_dir_path_obj, fixed = TRUE)) {
    celltypes <- c(
    '0' = "Oligodendrocyte",
    '1' = "Astrocyte", 
    '2' = "Ex. Neuron",
    '3' = "OPC Late",
    '4' = "Ex. Neuron-Signaling",#"Progenitor",
    '5' = "Ex. Neuron-Signaling",#"Progenitor",
    '6' = "Ex. Neuron",
    '7' = "In. Neuron",
    '8' = "Ex. Neuron-NEUROD6",
    '9' = "In. Neuron-TP53I11",
    '10' = "PVALB Interneuron",
    '11' = "Ex. Neuron-NEUROD6",
    '12' = "Ex. Neuron-NEUROD6",
    '13' = "Ex. Neuron",
    '14' = "Endothelia",
    '15' = "Ex. Neuron-NEUROD6",
    '16' = "Ex. Neuron-NEUROD6",
    '17' = "In. Neuron-TP53I11",
    '18' = "In. Neuron-BCL11B",
    '19' = "Ex. Neuron-NEUROD6",
    '20' = "Microglia",
    '21' = "Ex. Neuron-SATB2",
    '22' = "PVALB Interneuron",
    '23' = "Purkinje",
    '24' = "PVALB Interneuron",
    '25' = "Oligodendrocyte"
  )
} else if (grepl("gbm", parent_dir_path_obj, fixed = TRUE)) {
  
  if (grepl("oRG", filename, fixed = TRUE)) {
    celltypes <- c(
      '0' = 'Other',
      '1' = 'oRG'
    )
  } else {
    celltypes <- c(
      '0'  = "Mesenchymal", # "NPC (+TNC)",
      '1'  = "OPC Early", # "OPC Early",
      '2'  = "Astrocyte", # "Astrocyte", 
      '3'  = "Radial Glia 1", # "Astrocyte",
      '4'  = "Hypoxic", # "NPC (-TNC)",#"Neurons (hypoxic)",
      '5'  = "Microglia", # "Microglia", 
      '6'  = "Cycling 1", # "Cycling G1.S", 
      '7'  = "Mixed Cycling", # "Cycling (uncontrolled)", 
      '8'  = "Migrating Neuron", # "Neuron", 
      '9'  = "Pre-OPC", # "Pre-OPC", 
      '10' = "OPC Early Cycling", # "OPC G1.S",
      '11' = "Cycling 2", # "Cycling G2/M",
      '12' = "Oligodendrocyte", # "Oligodendrocyte", 
      '13' = "Immune", # "Immune",
      '14' = "Unknown", # "Unknown", # high MALAT1
      '15' = "OPC Late", # "OPC Late",
      '16' = "Radial Glia 2", # "Stressed Astrocyte", # with stress genes, less FAM107A
      '17' = "Endothelia", # "Endothelia",
      '18' = "Tumour Endothelia" # "Tumour Endothelia"
    )
  }
}

print(celltypes)
fccutoff <- 0.5
pcutoff <- 10e-5 
pmax <- 10e-32

cluster_markers$cluster_name <- ""
for (i in 1:length(celltypes)) {
  print(celltypes[i])
  cluster_markers[["cluster_name"]][which(cluster_markers[["cluster"]] == names(celltypes[i]))] <- celltypes[i]
}

if ("X" %in% colnames(cluster_markers)) {
  cluster_markers$gene <- sapply(strsplit(cluster_markers$X,split = "...", fixed = TRUE), function(x) x[1] )
  cluster_markers$X <- NULL
}

# set flag for p_value (for visualization)
cluster_markers$p_val_adj_limit <- ifelse(cluster_markers$p_val_adj <= pmax, pmax, cluster_markers$p_val_adj)
cluster_markers$flag <- cluster_markers$p_val_adj <= pmax #bonferroni default

# label genes by significance
cluster_markers <- cluster_markers %>% 
  mutate(fill_label = case_when(
            p_val_adj < pcutoff & avg_log2FC < -fccutoff ~ "p-value & Log2FC",
            p_val_adj < pcutoff & avg_log2FC >  fccutoff ~ "p-value & Log2FC",
            p_val_adj < pcutoff & (avg_log2FC < fccutoff & avg_log2FC > -fccutoff) ~ "p-value",
            p_val_adj >= pcutoff & avg_log2FC < -fccutoff ~ "Log2FC",
            p_val_adj >= pcutoff & avg_log2FC >  fccutoff ~ "Log2FC",
            .default = "NS"))
cluster_markers$fill_label <- factor(cluster_markers$fill_label, 
  levels=c("NS", "Log2FC", "p-value", "p-value & Log2FC"))
cluster_markers$p_val_log10_limit <- -log10(cluster_markers$p_val_adj_limit)

# join te_df
cluster_markers <- left_join(cluster_markers,te_df,by=c("gene"="gene_id"))

glimpse(cluster_markers)
write.csv(cluster_markers, file.path(figs_dir_path, paste0(cluster_group,"_tedeg.csv")))

cluster_markers <- filter(
  cluster_markers, 
  p_val_adj < pcutoff, 
  avg_log2FC < -fccutoff |  avg_log2FC > fccutoff, 
  type == "LINE" | type == "SINE" | type == "LTR" )

write.csv(cluster_markers, file.path(figs_dir_path, paste0(cluster_group,"_tedeg_filtered_10e-5p_0.5log2fc_types.csv")))

### End of Script
sessionInfo()
