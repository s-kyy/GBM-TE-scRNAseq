#!usr/bin/env Rscript

#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 1 filepath must be supplied: [xxx.rds]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1])){ 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    path_to_object <- dirname(obj_path)
    parent_dir_name <- basename(path_to_object)
  } else {
    stop("Filepath provided does not exist. Exiting...", call.=FALSE)
  }
}

#### Import Packages ####
set.seed(100)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(ggbreak)
library(scales) #plot axis manipulation

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr

set.seed(100)

#### Load Datasets ####
seurat.obj <- readRDS(obj_path) 
DefaultAssay(seurat.obj) <- "RNA"

if (is.factor(seurat.obj@meta.data$sample)) {
  sample_names <- sort(unique(droplevels(seurat.obj@meta.data$sample)))
} else {
  sample_names <- sort(unique(seurat.obj@meta.data$sample))
}
sample_names_indexed <- sample_names
names(sample_names_indexed) <- seq(1,length(sample_names))

#### Quality Control Figures ####
size <- 5
# Custom colours for up to 13 samples). 
samplePalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", 
                   "#ff716e", "#999999", "#0072B2", "#194c76", "#D55E00", 
                   "#3a4f41", "#6699cc", "#713e5a")

ifelse(!dir.exists(file.path(getwd(),parent_dir_name, "figs")),
        dir.create(file.path(getwd(),parent_dir_name, "figs"),recursive=T),
        "Directory Exists")

# Table of Cell Counts per Sample
df <- seurat.obj@meta.data %>% select(sample) %>% count(sample) 
write.csv(df, file= file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_samplecounts.csv")), row.names = F)

print("Making bar plot : cell count per sample")
p <- seurat.obj@meta.data %>%
    ggplot(aes(x=sample, fill=sample)) + 
    labs(fill="Sample ID") + labs(x = "", y = "Cell Count") +
    scale_fill_manual(values = samplePalette) + 
    geom_bar() +
    theme_classic() +
    # scale_y_break(c(30000, 70000), scales = 0.1, ticklabels=c(70000)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=10), legend.position = "none")
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_cellcounts.tiff")), 
       plot = p, units="in", width=size*1, height=size*0.6, dpi=300, compression = 'lzw')

print("Making density plot : Visuzlize number of UMIs or transcripts per cell")
p <- seurat.obj@meta.data %>% 
    ggplot(aes(color=sample, x=nUMI, fill=sample)) + 
    labs(x = "UMI per Transcript (log10)", y = "Log10 Cell Density") +
    labs(color="Sample ID") + scale_color_discrete(labels = sample_names) +
    # guides(fill=FALSE) + # Uncomment this line (and comment next line) if running <ggplot3.3.4
    guides(fill="none") + 
    geom_density(alpha = 0.2) + 
    scale_x_log10(labels = comma, limit=c(100, 100000)) + 
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = c(500,1000)) + geom_text(aes(x= 500, label="500\n", y = 1), angle=90, colour="black")
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_nUMI.tiff")), 
       plot = p, units="in", width=size*1.5, height=size, dpi=300, compression = 'lzw')

# Density Plot Visualize the distribution of genes detected per cell
p <- seurat.obj@meta.data %>% 
    ggplot(aes(color=sample, x=nGene, fill=sample)) + 
    labs(x = "Log10 Number of Genes Detected per Cell", y = "Log10 Cell Density") +
    labs(color="Sample ID") + scale_color_discrete(labels = sample_names) +
    # guides(fill=FALSE) + # Uncomment this line (and comment next line) if running <ggplot3.3.4
    guides(fill="none") + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10(labels = comma) + 
    geom_vline(xintercept = 300) +
    ggtitle("GrCh38+RT")
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_nGenes.tiff")), 
       plot = p, units="in", width=size*1.5, height=size, dpi=300, compression = 'lzw')

# Violin Plot: Visualize the distribution of genes detected per cell
p <- seurat.obj@meta.data %>% 
    ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
    labs(x = "Sample ID", y = "Log10 Number of Genes Detected per Cell") +
    labs(fill="Sample ID") + scale_fill_discrete(labels = sample_names) +
    geom_violin() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    ggtitle("GRCh38+RT")
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_nGenes_violin.tiff")), 
       plot = p, units="in", width=size*2, height=size, dpi=300, compression = 'lzw')

# Scatter Plot:  Visualize the correlation between genes detected and number of UMIs per cell
p <- seurat.obj@meta.data %>% 
    ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10(labels = comma) + 
    scale_y_log10() + 
    labs(x = "Log10 UMI per Transcript", y = "Log10 Genes Detected per Cell") +
    theme_classic() +
    geom_vline(xintercept = 500) + # Threshold for UMI per transcript. 
    geom_hline(yintercept = 250) + # Threshold for Genes detected per cell
    facet_wrap(~orig.ident, labeller = as_labeller(sample_names_indexed)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_nGenes_nUMI.tiff")), 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

# Density plot: Visualize the distribution of mitochondrial gene expression detected per cell
p <- seurat.obj@meta.data %>% 
    ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
    labs(x = "Mitochondrial Ratio", y = "Log10 Cell Density") +
    labs(color="Sample ID") + scale_color_discrete(labels = sample_names) +
    # guides(fill=FALSE) + # Uncomment this line (and comment next line) if running <ggplot3.3.4
    guides(fill="none") + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2) + 
    theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_mitoRatio.tiff")), 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

# Density plot Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
p <- seurat.obj@meta.data %>%
    ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
    labs(x = "Log10 Genes Detected Per UMI per Transcript", y = "Cell Count") +
    labs(color="Sample ID") + scale_color_discrete(labels = sample_names) +
    # guides(fill=FALSE) + # Uncomment this line (and comment next line) if running <ggplot3.3.4
    guides(fill="none") + 
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8) + 
    theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_novelGenes.tiff")), 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

#### End of Script ####
sessionInfo()
