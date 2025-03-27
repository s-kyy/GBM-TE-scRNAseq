#!usr/bin/env Rscript

library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("At least 2 filepaths must be supplied: [gte.rds] [te_df] [res] [figure_path]", call.=FALSE)
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
# library(reshape2)
library(dplyr)
library(viridis)
library(ggpubr) # stats in ggplot2

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(108)

#### ===================================================================== ####
#### Load Datasets ####
#### ===================================================================== ####
seurat_obj <- readRDS(obj_path) 
DefaultAssay(seurat_obj) <- "RNA"
te_df <- read.csv(te_path)
size    <- 5

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

glimpse(te_df)
glimpse(seurat_obj@meta.data)

#### ===================================================================== ####
#### Histogram of TE-types expressed across different celltypes ####
#### ===================================================================== ####

te_df$type <- factor(te_df$type, levels=c("LINE", "SINE", "LTR", "Other"))
p <- te_df  %>% ggplot(aes(x=num_celltypes, color=type, fill=type)) + 
  geom_histogram(binwidth=1, alpha=0.5) +
  labs(x="Number of Cell Types", y="Expressed TE subfamilies from Raw Counts",
    color="TE Type", fill="TE Type") +
  facet_grid(vars(type),scales="free")+
  scale_x_continuous(breaks=seq(0, max(te_df$num_celltypes)+1, by=1), 
    limits = c(0, max(te_df$num_celltypes)+1))+
  theme_classic()+
  theme(axis.text.x = element_text(size=12,color="black"), 
      axis.text.y = element_text(size=12,color="black"),
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      strip.text=element_text(size=12,face="bold"), 
      strip.background=element_rect(fill="grey",colour="grey"))
ggsave(file.path(figs_dir_path, paste0("hist_num_celltypes_splitby-TEtypes.tiff")),
    plot = p, units="in", width=size*1.3, height=size*0.9, dpi=300, compression = 'lzw')  


#### ===================================================================== ####
#### Percentage of counts per TE type per cell ####
#### ===================================================================== ####

# Barplot of TE-types per gsctypes
seurat_obj$LINE_perc <- PercentageFeatureSet(object=seurat_obj, features=te_df$gene_id[which((te_df$type == "LINE") )] )
seurat_obj$SINE_perc <- PercentageFeatureSet(object=seurat_obj, features=te_df$gene_id[which((te_df$type == "SINE") )] )
seurat_obj$LTR_perc <- PercentageFeatureSet(object=seurat_obj, features=te_df$gene_id[which((te_df$type == "LTR") )] )
seurat_obj$Other_perc <- PercentageFeatureSet(object=seurat_obj, features=te_df$gene_id[which((te_df$type == "Other") )] )

glimpse(seurat_obj@meta.data)
if ( "mitoRatio" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj@meta.data$mitoRatio <- NULL
}
write.csv(seurat_obj@meta.data, file.path(figs_dir_path, paste0("metadata_",filename,"_teRatios.csv")))
saveRDS(seurat_obj, file = file.path(parent_dir_path_obj, paste0(filename,"_teRatios.rds")))

#### ===================================================================== ####
#### Figures for Percentage of counts per TE type per cell ####
#### ===================================================================== ####

summary(seurat_obj$LINE_perc)
summary(seurat_obj$SINE_perc)
summary(seurat_obj$LTR_perc)
summary(seurat_obj$Other_perc)

makeUMAPPlot <- function(obj, feature) {
  DefaultAssay(obj) <- "RNA"
  p <- FeaturePlot(obj, features = feature)
  p <- p + theme(axis.line=element_blank(),
          axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  return(p)
}

p <- makeUMAPPlot(obj=seurat_obj, feature="LINE_perc")
ggsave(file.path(figs_dir_path, paste0("LINE_perc_FeaturePlot_nolabels.tiff")),
  plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj=seurat_obj, feature="SINE_perc")
ggsave(file.path(figs_dir_path, paste0("SINE_perc_FeaturePlot_nolabels.tiff")),
  plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj=seurat_obj, feature="LTR_perc")
ggsave(file.path(figs_dir_path, paste0("LTR_perc_FeaturePlot_nolabels.tiff")),
  plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj=seurat_obj, feature="Other_perc")
ggsave(file.path(figs_dir_path, paste0("Other_perc_FeaturePlot_nolabels.tiff")),
  plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')

makeUMAPPlot_cutoff <- function(obj, feature, cutoff=19.26) {
  DefaultAssay(obj) <- "RNA"
  p <- FeaturePlot(obj, features = feature,min.cutoff=cutoff)
  p <- p + theme(axis.line=element_blank(),
          axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  return(p)
}

p <- makeUMAPPlot_cutoff(obj=seurat_obj, feature="LINE_perc",cutoff=5)
ggsave(file.path(figs_dir_path, paste0("LINE_perc_FeaturePlot_nolabels_cutoff.tiff")),
  plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
p <- makeUMAPPlot_cutoff(obj=seurat_obj, feature="SINE_perc", cutoff=20)
ggsave(file.path(figs_dir_path, paste0("SINE_perc_FeaturePlot_nolabels_cutoff.tiff")),
  plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
# p <- makeUMAPPlot_cutoff(obj=seurat_obj, feature="LTR_perc", cutoff = 0)
# ggsave(file.path(figs_dir_path, paste0("LTR_perc_FeaturePlot_nolabels_cutoff.tiff")),
#   plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
# p <- makeUMAPPlot_cutoff(obj=seurat_obj, feature="Other_perc")
# ggsave(file.path(figs_dir_path, paste0("Other_perc_FeaturePlot_nolabels_cutoff.tiff")),
#   plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')


#### ===================================================================== ####
#### Violin Plots ####
#### ===================================================================== ####

if (grepl("gbm", parent_dir_path_obj, fixed = TRUE)) {
  violin_fill <- cluster_col_gbm_neftel
  violin_split <- cluster_col_gsc
  violin_colours <- c(
    "#FF6B6B", # "NEFTEL-AC",         pink-red
    "#A41917", # "NEFTEL-OPC",        dark red
    "#80DB84", # "NEFTEL-NPC1",       light green
    # "#daffde", # "NEFTEL-NPC2",       lighter green
    "#46a14a", # "NEFTEL-NPC1-2",     dark green
    "#1b761f", # "NEFTEL-NPC1-G1.S",  darker green
    "#99f7ff", # "NEFTEL-MES1",       lighter blue
    "#56B4E9", # "NEFTEL-MES2",       light blue
    "#2785ba", # "NEFTEL-MES1-2",     dark blue
    "#585191", # "NEFTEL-G1.S-G2.M",  darker blue
    "#A5B6AD" # "NA"  grey
  )
} else {
  violin_fill <- cluster_col
  violin_split <- cluster_col
}

for (te in c("LINE_perc", "SINE_perc", "LTR_perc", "Other_perc")) {
  te_label <- strsplit(te, "_", fixed=TRUE)[[1]][1]
  
  p <- seurat_obj@meta.data %>% ggplot(aes(x=!!sym(violin_split), y=!!sym(te), fill=!!sym(violin_fill))) + 
    geom_violin(trim=FALSE) + geom_boxplot(width=0.1, fill="white") + 
    labs(x="",y=paste("Percent of",te_label,"Counts (%)"), fill="GBM Subtype")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
    geom_hline(yintercept = median(seurat_obj@meta.data[[te]]), linetype = 2, color = "red")
  
  if (grepl("gbm", parent_dir_path_obj, fixed = TRUE)) {
    p <- p + scale_fill_manual(values=violin_colours)
  } else {
    p <- p + scale_fill_viridis_d()
  }
  
  ggsave(file.path(figs_dir_path, paste0(te,"_violin_teratio_nostats.tiff")),
    plot = p, units="in", width=size*1.5, height=size*0.7, dpi=300, compression = 'lzw')

  # p <- p + stat_compare_means(method = "kruskal.test", label.y = 85)+
  #   stat_compare_means(label = "p.format", method = "wilcox.test", ref.group = ".all.", 
  #     method.args=list(p.adjust.method = "BH"), hide.ns = TRUE, label.y = 80)  
  # ggsave(file.path(figs_dir_path, paste0(te,"_violin_teratio_stats.tiff")),
  #   plot = p, units="in", width=size*1.5, height=size*0.7, dpi=300, compression = 'lzw')

  p <- p + theme(legend.position = "none")
  ggsave(file.path(figs_dir_path, paste0(te,"_violin_teratio_noLegend.tiff")),
    plot = p, units="in", width=size*1.5, height=size*0.7, dpi=300, compression = 'lzw')
}

### End of Script
sessionInfo()
