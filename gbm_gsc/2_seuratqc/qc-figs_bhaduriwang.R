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
# library(ggbreak)
library(ggpmisc)
library(scales) #plot axis manipulation

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr

set.seed(100)

#### Load Datasets ####
seurat.obj <- readRDS(obj_path) 
DefaultAssay(seurat.obj) <- "RNA"

# Factor sample names
if (is.factor(seurat.obj@meta.data$orig.ident)) {
  seurat.obj@meta.data$orig.ident <- droplevels(seurat.obj@meta.data$orig.ident)
}
if (is.factor(seurat.obj@meta.data$sample)) {
  sample_names <- sort(unique(droplevels(seurat.obj@meta.data$sample)))
} else {
  sample_names <- sort(unique(seurat.obj@meta.data$sample))
}
seurat.obj@meta.data$sample <- factor(
  seurat.obj@meta.data$sample, 
  levels = sample_names
)
sample_names_indexed <- sample_names
names(sample_names_indexed) <- seq(1,length(sample_names))

#### Quality Control Figures ####
size <- 5
# Custom colours for up to 13 samples). 
# sample_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", 
#                    "#ff716e", "#999999", "#0072B2", "#194c76", "#D55E00", 
#                    "#3a4f41", "#6699cc", "#713e5a")
sample_palette<- c(
  "#E69F00", "#56B4E9", "#009E73", 
  "#F0E442", "#CC79A7", "#ff716e",
  "#999999", "#0072B2", "#194c76")

ifelse(!dir.exists(file.path(getwd(),parent_dir_name, "figs")),
        dir.create(file.path(getwd(),parent_dir_name, "figs"),recursive=T),
        "Directory Exists")

print("Bar plot Cell Count per Sample")
p <- seurat.obj@meta.data %>%
    ggplot(aes(x=sample, fill=sample)) + 
    labs(x = "", y = "Cell Count") +
    scale_fill_manual(values = sample_palette) + 
    labs(fill="Sample ID") + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=10), legend.position = "none")
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_cellcounts.tiff")), 
       plot = p, units="in", width=size*1, height=size*0.8, dpi=300, compression = 'lzw')

print("Density Plot: Visualize the number UMIs/transcripts per cell")
p <- seurat.obj@meta.data %>% 
    ggplot(aes(x=nUMI, color=sample, fill=sample)) + #y=..count.. / ..scaled are other options fro geom_density()
    labs(x = "Log10 Number of UMI/Transcripts Detected", y = "Cell Density") +
    scale_color_manual(values = sample_palette) +
    scale_fill_manual(values = sample_palette) +
    labs(fill="Sample ID", color="Sample ID") + 
    geom_density(alpha = 0.2) + 
    scale_x_log10(labels = comma) + #, limits=c(0,100000)
    theme_classic() +
    geom_vline(xintercept = c(1000), colour = "red") #+
    # geom_text(aes(x= 1000, label="1000\n", y = 1), angle=90, colour="red")
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_nUMI_log10x.tiff")), 
       plot = p, units="in", width=size*1.1, height=size*0.8, dpi=300, compression = 'lzw')

print("Density Plot: Visualize the number UMIs/transcripts per cell")
p <- seurat.obj@meta.data %>% 
    ggplot(aes(x=nUMI, color=sample, fill=sample)) + #y=..count.. / ..scaled are other options fro geom_density()
    labs(x = "Number of UMI/Transcripts Detected", y = "Cell Density") +
    scale_color_manual(values = sample_palette) +
    scale_fill_manual(values = sample_palette) +
    labs(fill="Sample ID", color="Sample ID") + 
    geom_density(alpha = 0.2) + 
    # scale_x_log10(labels = comma) + #, limits=c(0,100000)
    theme_classic() #+
    # geom_vline(xintercept = c(1000), colour = "red") +
    # geom_text(aes(x= 1000, label="1000\n", y = 1), angle=90, colour="red")
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_nUMI.tiff")), 
       plot = p, units="in", width=size*1.1, height=size*0.8, dpi=300, compression = 'lzw')

print("Density Plot Visualize the distribution of genes detected per cell")
p <- seurat.obj@meta.data %>% 
    ggplot(aes(x=nGene, color=sample, fill=sample)) + 
    labs(x = "Log10 Number of Genes Detected", y = "Cell Density") +
    scale_color_manual(values = sample_palette) +
    scale_fill_manual(values = sample_palette) +
    labs(fill="Sample ID", color="Sample ID") + 
    scale_x_log10(labels = comma) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    geom_vline(xintercept = 300, colour="red") 
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_nGenes_log10x.tiff")), 
       plot = p, units="in", width=size*1.1, height=size*0.8, dpi=300, compression = 'lzw')

print("Density Plot Visualize the distribution of genes detected per cell")
p <- seurat.obj@meta.data %>% 
    ggplot(aes(x=nGene, color=sample, fill=sample)) + 
    labs(x = "Number of Genes Detected", y = "Cell Density") +
    scale_color_manual(values = sample_palette) +
    scale_fill_manual(values = sample_palette) +
    labs(fill="Sample ID", color="Sample ID") + 
    # scale_x_log10(labels = comma) + 
    geom_density(alpha = 0.2) + 
    theme_classic() #+
    # geom_vline(xintercept = 300, colour="red") +
    # geom_text(aes(x= 300, label="300\n", y = 1), angle=90, colour="red")
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_nGenes.tiff")), 
       plot = p, units="in", width=size*1.1, height=size*0.8, dpi=300, compression = 'lzw')
my.formula <- y ~ x 
print("Scatter Plot:  Visualize the correlation between genes detected and number of UMIs per cell")
p <- seurat.obj@meta.data %>% 
    ggplot(aes(x=nUMI, y=nGene, color=mitoRatio), group=sample) + 
    labs(x = "Log10 Number of UMI/Transcripts Detected", y = "Log10 Number of Genes Detected") +
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    # stat_smooth(method=lm) +
    # stat_poly_line() +
    # stat_poly_eq() + 
    geom_smooth(method = "lm", se=FALSE, formula = my.formula, color="blue") +
    stat_poly_eq(geom = "text_npc",
               formula = my.formula, parse = TRUE,
              #  aes(label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),
               aes(label =  paste(..rr.label..)),
               color="blue", size=3) +
    scale_x_log10(labels = comma) + 
    scale_y_log10(labels = comma) + 
    labs(color="Mito. Ratio") +
    theme_classic() +
    geom_vline(xintercept = 1000, linetype = "dashed", colour = "red") +
    geom_hline(yintercept = 300, linetype = "dashed", colour = "red") +
    facet_wrap(~sample) + 
    # facet_wrap(~sample, labeller = as_labeller(sample_names_indexed)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
      axis.title.y = element_text(size=11, hjust=0.5), 
      strip.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black"))
    # theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_nGenes_nUMI_mitoRatio_log10x.tiff")), 
       plot = p, units="in", width=size*1.2, height=size*1.2, dpi=300, compression = 'lzw')

print("Density plot: Visualize the distribution of mitochondrial gene expression detected per cell")
p <- seurat.obj@meta.data %>% 
    ggplot(aes(x=mitoRatio, color=sample, fill=sample)) + 
    labs(x = "Mitochondrial Ratio", y = "Cell Density") +
    scale_color_manual(values = sample_palette) +
    scale_fill_manual(values = sample_palette) +
    labs(fill="Sample ID", color="Sample ID") + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    geom_vline(xintercept = 0.2, colour = "red") #+ 
    # geom_text(aes(x= 0.2, label="0.2\n", y = 10), angle=90, colour="red")
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_mitoRatio.tiff")), 
       plot = p, units="in", width=size*1, height=size*0.8, dpi=300, compression = 'lzw')

print("Density plot Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI")
p <- seurat.obj@meta.data %>%
    ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
    labs(x = "Log10 Genes Detected Per UMI per Transcript", y = "Cell Density") +
    scale_color_manual(values = sample_palette) +
    scale_fill_manual(values = sample_palette) +
    labs(fill="Sample ID", color="Sample ID") + 
    geom_density(alpha = 0.2) +
    theme_classic() +
    xlim(0.75, 0.96)+
    geom_vline(xintercept = 0.8, colour = "red") 
ggsave(file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_novelGenes.tiff")), 
       plot = p, units="in", width=size*1, height=size*0.8, dpi=300, compression = 'lzw')
rm("p")

# Table of Cell Counts per Sample
df <- seurat.obj@meta.data %>% 
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
obj.list <- SplitObject(seurat.obj, split.by="sample")
rm("seurat.obj")
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
write.csv(df, file= file.path(getwd(), parent_dir_name, "figs", paste0(filename,"_samplecounts.csv")), row.names = F)

#### End of Script ####
sessionInfo()
