#!usr/bin/env Rscript

# .libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
# .libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 2 filepath must be supplied: [markers.csv] [fig_dir_name]", call.=FALSE)
} else if (length(args)>=1) {
  
  # verify filepaths
  if (file.exists(args[1]) ) { 
    gsea_path <- args[1]
    filename <- basename(path_ext_remove(gsea_path))
    parent_dir_path_obj <- dirname(gsea_path)
  } else {
    stop("Filepaths provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>1) {
    figs_dir_name <- args[2]
  } else {
    figs_dir_name <- "figs_gsea"
  }

} else {
    stop("Error in calling R or filepaths. Exiting...", call.=FALSE)
}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
set.seed(108)

library(ggplot2)
library(tidyverse) 
library(dplyr)
library(viridis)

set.seed(108)
options(warn=1) #print warning messages as they occur

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
gsea_df <- read.csv(gsea_path)
size <- 5

figs_dir_path <- file.path(parent_dir_path_obj,figs_dir_name)

ifelse(!dir.exists(file.path(figs_dir_path)),
        dir.create(file.path(figs_dir_path),recursive=T),
        "Directory Exists")

#### =========================================== ####
#### Generate GSEA Figures ####
#### =========================================== ####

glimpse(gsea_df)

# horizontal bar plot with nlog10(FDR) x Geneset labeled by Size on secondary axis
# Columns split by celltype # Rows split by geneset category

wrap_label <- function(x) str_wrap(str_replace_all(x, "-", " "), width = 3)

for (i in sort(unique(gsea_df$Cluster))) {
  
  print(paste("Making plot",i))
  
  # Filter dataframe
  temp <- gsea_df %>% dplyr::filter(Cluster == i) %>% 
    mutate(across(c(NES), \(x) round(x,2))) %>% 
    mutate(Geneset = reorder(Geneset, nlog10p, decreasing = FALSE))
  
  # Store number of terms to dynamically adjust figure size
  num_terms <- length(temp[,1])
  if (num_terms <= 5) num_terms = 6 
  
  # Generate bar plot
  p <- temp %>% 
    ggplot(aes(x=nlog10p, y=Geneset,fill=nlog10p)) + 
    geom_bar(stat="identity") +
    labs(x="", y="", fill="-log10(FDR q-value)")+
    facet_grid(Collection ~ . ,scales="free",space = "free",labeller = as_labeller(wrap_label))+
    geom_text(aes(label = NES), colour = "white",hjust=1.5, fontface = "bold", size=4)+
    theme_minimal()+
    theme(axis.text = element_text(size = 13),
      strip.text=element_text(size=10,face="bold"), 
      strip.background=element_rect(fill="grey",colour="grey"),
      legend.position = "bottom")
  ggsave(file.path(figs_dir_path, paste0(filename,"_",i,"bar.tiff")),
      plot = p, units="in", width=size*1.7, height=size*num_terms*(0.8/7.0), dpi=300, compression = 'lzw')
}

if ("significance" %in% names(gsea_df)) {
  wrap_label <- function(x) str_wrap(gsub('[0-9]+', '', x), width = 5)
  num_terms <- length(unique(gsea_df[,1]))
  if (num_terms <= 5) num_terms = 6 

  p <- gsea_df %>% ggplot(aes(x=factor(Cluster), y=Geneset,size=NES,colour=factor(significance, levels=c("p<0.05","p<0.01","p<1.0E-6")))) + 
    # geom_bar(stat="identity") +
    geom_point() +
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    labs(x="Cluster", y="", colour="FDR q-value")+
    facet_grid(Collection ~ . ,scales="free",space = "free",labeller = as_labeller(wrap_label))+
    geom_text(aes(label = Size), colour = "black", vjust=-1.25, fontface = "bold", size=2.5)+
    theme_minimal() +
    theme(axis.text = element_text(size = 13),
      strip.text=element_text(size=10,face="bold"), 
      strip.background=element_rect(fill="grey",colour="grey"),
      legend.position = "bottom")
  ggsave(file.path(figs_dir_path, paste0(filename,"_full_bar.tiff")),
      plot = p, units="in", width=size*1.6, height=size*num_terms*0.1, dpi=300, compression = 'lzw')
}

## End of script
sessionInfo()