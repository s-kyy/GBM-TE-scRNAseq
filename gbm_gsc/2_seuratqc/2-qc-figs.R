#!usr/bin/env Rscript

#### Import Packages ####
set.seed(34)

library(Seurat)
library(Matrix)
library(ggplot2)
library(genefilter)
library(tidyverse)
library(scales)
library(AnnotationHub)

set.seed(34)

#### Load Datasets ####
sampleNames <- read.csv("./samples.csv", header = FALSE)  # USER ADJUST
sampleNames <- sampleNames[,1]
sampleNames_named <- sampleNames
names(sampleNames_named) <- seq(1,length(sampleNames))
gte <- readRDS("temp/gte.rds") # USER ADJUST
ge <- readRDS("temp/ge.rds") # USER ADJUST
DefaultAssay(gte) <- "RNA"
DefaultAssay(ge) <- "RNA"

#### Quality Control Figures ####
size <- 7

p <- gte@meta.data %>%
    ggplot(aes(x=orig.ident, fill=orig.ident)) + 
    labs(fill="samples") + labs(x = "samples", y = "cell count") +
    scale_x_discrete(labels = sampleNames) +
    scale_fill_discrete(labels = sampleNames) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GRCh38+RT") + theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave("temp/gte_cellcounts.tiff", 
       plot = p, units="in", width=size*1.5, height=size, dpi=300, compression = 'lzw')


p <- ge@meta.data %>%
    ggplot(aes(x=orig.ident, fill=orig.ident)) + 
    labs(fill="samples") + labs(x = "samples", y = "cell count") +
    scale_x_discrete(labels = sampleNames) +
    scale_fill_discrete(labels = sampleNames) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GRCh38") + theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave("temp/ge_cellcounts.tiff", 
       plot = p, units="in", width=size*1.5, height=size, dpi=300, compression = 'lzw')

# Visualize the number UMIs/transcripts per cell
p <- gte@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nUMI, fill=orig.ident)) + 
    labs(x = "UMI per Transcript (log10)", y = "log 10 Cell Density") +
    labs(color="samples") + scale_color_discrete(labels = sampleNames) +
    guides(fill=FALSE) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10(labels = comma, limit=c(100, 100000)) + 
    theme_classic() +
    ggtitle("GRCh38+RT") + theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = c(500,1000)) + geom_text(aes(x= 500, label="500\n", y = 1), angle=90, colour="black")
ggsave("temp/gte_nUMI.tiff", 
       plot = p, units="in", width=size*1.5, height=size, dpi=300, compression = 'lzw')

p <- ge@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nUMI, fill=orig.ident)) + 
    labs(x = "UMI per Transcript (log10)", y = "log 10 Cell Density") +
    labs(color="samples") + scale_color_discrete(labels = sampleNames) +
    guides(fill=FALSE) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10(labels = comma, limit=c(100, 100000)) + 
    theme_classic() +
    ggtitle("GRCh38") + theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = c(500,1000)) + geom_text(aes(x= 500, label="500\n", y = 1), angle=90, colour="black")
ggsave("temp/ge_nUMI.tiff", 
       plot = p, units="in", width=size*1.5, height=size, dpi=300, compression = 'lzw')

# Visualize the distribution of genes detected per cell
p <- gte@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nGene, fill=orig.ident)) + 
    labs(x = "log 10 Number of Genes Detected per Cell", y = "log 10 Cell Density") +
    labs(color="samples") + scale_color_discrete(labels = sampleNames) +
    guides(fill=FALSE) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10(labels = comma) + 
    geom_vline(xintercept = 300) +
    ggtitle("GrCh38+RT")
ggsave("temp/gte_nGenes.tiff", 
       plot = p, units="in", width=size*1.5, height=size, dpi=300, compression = 'lzw')

p <- ge@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nGene, fill=orig.ident)) + 
    labs(x = "log 10 Number of Genes Detected per Cell", y = "log 10 Cell Density") +
    labs(color="samples") + scale_color_discrete(labels = sampleNames) +
    guides(fill=FALSE) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10(labels = comma) + 
    geom_vline(xintercept = 300) +
    ggtitle("GrCh38")
ggsave("temp/ge_nGenes.tiff", 
       plot = p, units="in", width=size*1.5, height=size, dpi=300, compression = 'lzw')

# Visualize the distribution of genes detected per cell via violin plot
p <- gte@meta.data %>% 
    ggplot(aes(x=orig.ident, y=log10(nGene), fill=orig.ident)) + 
    labs(x = "samples", y = "log 10 Number of Genes Detected per Cell") +
    labs(fill="samples") + scale_fill_discrete(labels = sampleNames) +
    geom_violin() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    ggtitle("GRCh38+RT")
ggsave("temp/gte_nGenes_violin.tiff", 
       plot = p, units="in", width=size*2, height=size, dpi=300, compression = 'lzw')

p <- ge@meta.data %>% 
    ggplot(aes(x=orig.ident, y=log10(nGene), fill=orig.ident)) + 
    labs(x = "samples", y = "log 10 Number of Genes Detected per Cell") +
    labs(fill="samples") + scale_fill_discrete(labels = sampleNames) +
    geom_violin() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    ggtitle("GRCh38")
ggsave("temp/ge_nGenes_violin.tiff", 
       plot = p, units="in", width=size*2, height=size, dpi=300, compression = 'lzw')

# Visualize the correlation between genes detected and number of UMIs 
p <- gte@meta.data %>% 
    ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10(labels = comma) + 
    scale_y_log10() + 
    labs(x = "log 10 UMI per Transcript", y = "log 10 Genes Detected per Cell") +
    theme_classic() +
    geom_vline(xintercept = 500) + # Threshold for UMI per transcript. 
    geom_hline(yintercept = 250) + # Threshold for Genes detected per cell
    facet_wrap(~orig.ident, labeller = as_labeller(sampleNames_named)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GRCh38+RT") + theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave("temp/gte_nGenes_nUMI.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

p <- ge@meta.data %>% 
    ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10(labels = comma) + 
    scale_y_log10() + 
    labs(x = "log 10 UMI per Transcript", y = "log 10 Genes Detected per Cell") +
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~orig.ident, labeller = as_labeller(sampleNames_named)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("GRCh38") + theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave("temp/ge_nGenes_nUMI.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

# Visualize the distribution of mitochondrial gene expression detected per cell
p <- gte@meta.data %>% 
    ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
    labs(x = "Mitochondrial Ratio", y = "log 10 Cell Density") +
    labs(color="samples") + scale_color_discrete(labels = sampleNames) +
    guides(fill=FALSE) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2) + 
    ggtitle("GRCh38+RT") + theme(plot.title = element_text(hjust = 0.5))
ggsave("temp/gte_mitoRatio.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

p <- ge@meta.data %>% 
    ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
    labs(x = "Mitochondrial Ratio", y = "log 10 Cell Density") +
    labs(color="samples") + scale_color_discrete(labels = sampleNames) +
    guides(fill=FALSE) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2) + 
    ggtitle("GRCh38") + theme(plot.title = element_text(hjust = 0.5))
ggsave("temp/ge_mitoRatio.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
p <- gte@meta.data %>%
    ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
    labs(x = "log 10 Genes Detected Per UMI per Transcript", y = "Cell Count") +
    labs(color="samples") + scale_color_discrete(labels = sampleNames) +
    guides(fill=FALSE) + 
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8) + 
    ggtitle("GRCh38+RT") + theme(plot.title = element_text(hjust = 0.5))
ggsave("temp/gte_novelGenes.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

p <- ge@meta.data %>%
    ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
    labs(x = "log 10 Genes Detected Per UMI per Transcript", y = "Cell Count") +
    labs(color="samples") + scale_color_discrete(labels = sampleNames) +
    guides(fill=FALSE) + 
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8) + 
    ggtitle("GRCh38") + theme(plot.title = element_text(hjust = 0.5))
ggsave("temp/ge_novelGenes.tiff", 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

#### End of Script ####
sessionInfo()
