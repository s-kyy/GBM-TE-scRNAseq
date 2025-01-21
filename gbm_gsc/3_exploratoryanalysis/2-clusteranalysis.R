#!usr/bin/env Rscript

# preliminary analysis of celltypes and marker volcano plots

library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("At least 3 filepaths must be supplied: [xxx.rds] [xxx_markers_0.3.csv] [xxx_markers_0.4.csv] [figs_dir_name] [res1] [res2]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2]) && file.exists(args[3])) { 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_name_obj <- basename(parent_dir_path_obj)
    marker_path3 <- args[2]
    marker_path4 <-args[3]
#     parent_dir_path_marker <- dirname(marker_path)
#     parent_dir_name_marker <- basename(parent_dir_path_marker)
  } else {
    stop("Filepaths provided do not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>3) {
    figs_dir_name <- args[4]
  } else {
    figs_dir_name <- "figs_clusteranalysis"
  }

  if (length(args)>4) {
    cluster_res03 <- args[5]
  } else {
    cluster_res03 <- "integrated_snn_res.0.3"
  }

  if (length(args)>5) {
    cluster_res04 <- args[6]
  } else {
    cluster_res04 <- "integrated_snn_res.0.4"
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
library(reshape2)
library(viridis)
library(ggrepel)

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(108)

#### ===================================================================== ####
#### Load Datasets ####
#### ===================================================================== ####
seurat.obj <- readRDS(obj_path) 
cluster_markers3 <- read.csv(marker_path3)
cluster_markers4 <- read.csv(marker_path4)
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

seurat.obj@meta.data[[cluster_res03]] <- factor(seurat.obj@meta.data[[cluster_res03]], levels= sort(unique(seurat.obj@meta.data[[cluster_res03]])))
seurat.obj@meta.data[[cluster_res04]] <- factor(seurat.obj@meta.data[[cluster_res04]], levels= sort(unique(seurat.obj@meta.data[[cluster_res04]])))

figs_dir_path <- file.path(parent_dir_path_obj, figs_dir_name)

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

print(seurat.obj)

#### ===================================================================== ####
#### UMAP plots by cluster ####
#### ===================================================================== ####
DefaultAssay(seurat.obj) <-"RNA"

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res03) +
              scale_colour_viridis_d()
ggsave(file.path(figs_dir_path, paste0(cluster_res03,"_UMAP_nolabels.tiff")),
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res03, 
              label = TRUE, 
              label.size = 3.5) +
              scale_colour_viridis_d()
ggsave(file.path(figs_dir_path, paste0(cluster_res03,"_UMAP.tiff")),
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res04) +
              scale_colour_viridis_d()
ggsave(file.path(figs_dir_path, paste0(cluster_res04,"_UMAP_nolabel.tiff")),
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res04, 
              label = TRUE, 
              label.size = 3.5) +
              scale_colour_viridis_d()
ggsave(file.path(figs_dir_path, paste0(cluster_res04,"_UMAP.tiff")),
       plot = p, units="in", width=size*1.1, height=size*1, dpi=300, compression = 'lzw')

print("UMAPs by cluster exported")

#### ===================================================================== ####
#### UMAP plots by sample ####
#### ===================================================================== ####

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res03,
              split.by = "sample") +
              scale_colour_viridis_d()
ggsave(file.path(figs_dir_path, paste0("sample_split_", cluster_res03,"_UMAP.tiff")),
       plot = p, units="in", width=size*6, height=size*0.8, dpi=300, compression = 'lzw')

p <- DimPlot(seurat.obj, reduction = "umap", 
              group.by = cluster_res04,
              split.by = "sample") +
              scale_colour_viridis_d()
ggsave(file.path(figs_dir_path, paste0("sample_split_", cluster_res04,"_UMAP.tiff")),
       plot = p, units="in", width=size*6, height=size*0.8, dpi=300, compression = 'lzw')

p <- DimPlot(seurat.obj, reduction = "umap", 
      group.by = "sample", cols = sample_palette) 
ggsave(file.path(figs_dir_path, paste0("sample_custom_cols","_UMAP.tiff")),
      plot = p, units="in", width=size*1.3, height=size*1, dpi=300, compression = 'lzw')
print("UMAPs by sample_orig exported")

if (grepl("gbm", parent_dir_name_obj, fixed = TRUE)) {
       
  p <- DimPlot(seurat.obj, reduction = "umap", 
        group.by = "sample_orig") 
  ggsave(file.path(figs_dir_path, paste0("sample_orig","_UMAP.tiff")),
        plot = p, units="in", width=size*1.3, height=size*1, dpi=300, compression = 'lzw')
  print("UMAPs by sample_orig exported")

  p <- DimPlot(seurat.obj, reduction = "umap", 
                group.by = cluster_res03,
                split.by = "sample_orig") +
              scale_colour_viridis_d()
  ggsave(file.path(figs_dir_path, paste0("sample_orig_split_", cluster_res03,"_UMAP.tiff")),
         plot = p, units="in", width=size*7, height=size*0.8, dpi=300, compression = 'lzw')
  
  p <- DimPlot(seurat.obj, reduction = "umap", 
                group.by = cluster_res04,
                split.by = "sample_orig") +
              scale_colour_viridis_d()
  ggsave(file.path(figs_dir_path, paste0("sample_orig_split_", cluster_res04,"_UMAP.tiff")),
         plot = p, units="in", width=size*7, height=size*0.8, dpi=300, compression = 'lzw')
  
  print("UMAPs split by sample exported")
}

#### ===================================================================== ####
#### UMAP plots by female sample ####
#### ===================================================================== ####
p <- DimPlot(seurat.obj, 
              reduction = "umap", 
              cells.highlight = WhichCells(object = seurat.obj, 
              expression = (sample == female_samples) ))
ggsave(file.path(figs_dir_path, paste0("female_samples","_UMAP.tiff")),
      plot = p, units="in", width=size*1.4, height=size*1, dpi=300, compression = 'lzw')

print("UMAPs by female_sample exported")

#### ===================================================================== ####
#### Ridge Plots ####
#### ===================================================================== ####
DefaultAssay(seurat.obj) <- "integrated"

Idents(seurat.obj) <- cluster_res03
p <- RidgePlot(seurat.obj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
ggsave(file.path(figs_dir_path, paste0( cluster_res03,"_cellcycleRidgePlot_integrated.tiff")), 
      plot = p, units="in", width=size*1.2, height=size*1.3, dpi=300, compression = 'lzw')
Idents(seurat.obj) <- cluster_res04
p <- RidgePlot(seurat.obj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
ggsave(file.path(figs_dir_path, paste0( cluster_res04,"_cellcycleRidgePlot_integrated.tiff")), 
      plot = p, units="in", width=size*1.2, height=size*1.3, dpi=300, compression = 'lzw')

#### ===================================================================== ####
#### Differential Gene Expression Figures ####
#### ===================================================================== ####

#### MA plot (instead of Volcano plot) 
# help visualize highly expressed genes per cluster that have 
# a p-value too small for R to export. 

## 1-Calculate mean expression of marker gene per cluster and add to table
## Save table
meanExpressionPerCluster <- function(seurat.object, cluster_df, idents) {

  DefaultAssay(seurat.object) <- "RNA"
  seurat.object <- NormalizeData(seurat.object) #default: LogNormalization

  # calculate mean expression of genes per cluster
  avg_exp_df <- AverageExpression(
    seurat.obj,
    assays = "RNA",
    slot = "data",
    features = unique(cluster_df$gene),
    group.by = idents
  )

  avg_exp_df <- as.data.frame(avg_exp_df$RNA)
  avg_exp_df <- rownames_to_column(avg_exp_df, var = "gene")
  avg_exp_df <- melt(avg_exp_df, value.name="avg.exp", variable.name="cluster")
  avg_exp_df$cluster <- as.numeric(as.character(avg_exp_df$cluster))
  cluster_df <- left_join(x = cluster_df, y = avg_exp_df, 
                          by = c("cluster", "gene"))
  cluster_df$p_val.bonf <- p.adjust(cluster_df$p_val, 
                                          method = "bonferroni")
  cluster_df <- cluster_df %>% 
    mutate(fill_label = case_when(
              p_val.bonf < 10e-20 ~ "p.adj < 10e-20",
              p_val.bonf < 0.01 ~ "p.adj < 0.01",
              p_val.bonf < 0.05 ~ "p.adj < 0.05",
              p_val.bonf >= 0.05 ~ "p.adj >= 0.05"))
  cluster_df$avg.exp.limit <- ifelse(cluster_df$avg.exp > 100, 100, cluster_df$avg.exp)
  cluster_df$avg.exp.flag <- cluster_df$avg.exp > 100
  print(head(cluster_df$fill_label))
  print(head(cluster_df$gene))
  cluster_df$fill_label <- factor(cluster_df$fill_label, levels=c("p.adj < 10e-20", "p.adj < 0.01", "p.adj < 0.05", "p.adj >= 0.05"))
  return(cluster_df)  
}

## 2- generate MA plot with Mean Normalized Expression as x -log2(FC) as y
## impossible p-values = p < 10e-20
## below 0.01 / 0.05
## above 0.05
MAPlot <- function(cluster_df,cols,resolution, min_exp_threshold) {
  cluster_id <- unique(unlist(cluster_df$cluster))
  print(cluster_id)
  for (i in 1:length(cluster_id)){
    cluster_df_subset <- subset(cluster_df, cluster == cluster_id[i])
    p <-  ggplot() + #, group=cluster) +
          labs(x = "Mean Normalized Expression", y = "Log2(Fold-Change)") +
          geom_point(data = subset(cluster_df_subset, !avg.exp.flag),
                aes(x=avg.exp, y = avg_log2FC, color=fill_label), size=2)+
          geom_point(data = subset(cluster_df_subset, avg.exp.flag),
                aes(x=avg.exp.limit, y = avg_log2FC, color = "indianred"))+
          scale_color_manual(values=c("seagreen2", "turquoise", "royalblue2", "grey", "indianred")) +
          coord_cartesian(xlim = c(-0.5, 125), clip="off") +
          theme_minimal() + 
          theme(strip.background = element_blank(),
            axis.line.y = element_line(colour="grey50"),
            legend.title = element_blank(),
            legend.position.inside = c(.95, .95),
            legend.justification = c("right", "bottom"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6),
            legend.background = element_rect(fill = "white", colour = "white"))
    print(paste("Saving figure of cluster",cluster_id[i]))
    ggsave(file.path(figs_dir_path, paste0(resolution,"_",cluster_id[i],"_MAPlot.tiff")),
      plot = p, units="in", width=size*1, height=size*0.8, dpi=300, compression = 'lzw')

    p <- p + # theme(plot.margin=unit(c(1, 3, 1, 1), "cm")) +
        geom_text_repel(
          data = head(subset(cluster_df_subset, 
                        # (avg_log2FC > 1 | avg_log2FC < -1) & 
                        p_val.bonf < 0.01 & 
                        avg.exp > min_exp_threshold & avg.exp < 100), 20),
          mapping = aes(x = avg.exp, y = avg_log2FC, label = gene),
          max.overlaps = Inf, seed = 34, force = 2,
          nudge_x = 5, min.segment.length = 0.1, size = 2,
          box.padding = 0.5, xlim = c(-0.5, 125)) + 
        geom_text_repel(
          data = subset(cluster_df_subset, 
                        # (avg_log2FC > 1 | avg_log2FC < -1) & 
                        p_val.bonf < 0.01 & avg.exp >= 100),
          mapping = aes(x = avg.exp.limit, y = avg_log2FC, label = gene),
          max.overlaps = Inf, seed = 34, force = 2, direction = "y",
          nudge_x = 10, min.segment.length = 0.1, size = 2,
          box.padding = 0.5, xlim = c(-0.5, 125))  

    print(paste("Saving figure of cluster",cluster_id[i], "with labels"))
    
    ggsave(file.path(figs_dir_path, paste0(resolution,"_",cluster_id[i],"_MAPlot_labeled.tiff")),
      plot = p, units="in", width=size*1, height=size*0.8, dpi=300, compression = 'lzw')
  }
}

DefaultAssay(seurat.obj) <- "RNA"

if ("X" %in% colnames(cluster_markers3) &&  "X" %in% colnames(cluster_markers4)) {
  cluster_markers3$gene <- sapply(strsplit(
                            cluster_markers3$X,split = "...", fixed = TRUE), function(x) x[1] )
  cluster_markers4$gene <- sapply(strsplit(
                            cluster_markers4$X,split = "...", fixed = TRUE), function(x) x[1] )
  
  cluster_markers3$X <- NULL
  cluster_markers4$X <- NULL
}

cluster_markers3 <- meanExpressionPerCluster(
  seurat.object = seurat.obj,
  cluster_df = cluster_markers3,
  idents = cluster_res03
)
MAPlot(cluster_markers3, p_val_cols, cluster_res03, avg_exp_threshold)

cluster_markers4 <- meanExpressionPerCluster(
  seurat.object = seurat.obj,
  cluster_df = cluster_markers4,
  idents = cluster_res04
)
MAPlot(cluster_markers4, p_val_cols, cluster_res04,avg_exp_threshold)

#### ===================================================================== ####
#### Create Heatmap ####
#### ===================================================================== ####

# #make heatmap
# Heatmapplot <- function(seurat.object, feature_list, idents) {
#   DefaultAssay(seurat.object) <- "RNA"
#   p <- DoHeatmap (seurat.object,
#                 slot = "data",
#                 group.by = idents,
#                 features = feature_list)
#                 # features = factor(unique(cluster_df$gene), 
#                 #                   levels = unique(cluster_df$gene)) ) 
#   return(p)
# }

# Filter out genes p-val.bonf >= 0.01 
print(dim(cluster_markers3))
cluster_markers3 <- cluster_markers3 %>% 
  # filter((avg_log2FC > 1 | avg_log2FC < -1) & p_val.bonf < 0.01 )
  filter(p_val.bonf < 0.01 )
print(dim(cluster_markers3))
# write.csv(cluster_markers3, file.path(figs_dir_path, paste0(cluster_res03,"_log2FC1_padj0.01_DEG.csv")), row.names=FALSE)
write.csv(cluster_markers3, file.path(figs_dir_path, paste0(cluster_res03,"_padj0.01_DEG.csv")), row.names=FALSE)

print(dim(cluster_markers4))
cluster_markers4 <- cluster_markers4 %>% 
  # filter((avg_log2FC > 1 | avg_log2FC < -1) & p_val.bonf < 0.01)
  filter(p_val.bonf < 0.01)
print(dim(cluster_markers4))
# write.csv(cluster_markers4, file.path(figs_dir_path, paste0(cluster_res04,"_log2FC1_padj0.01_DEG.csv")), row.names=FALSE)
write.csv(cluster_markers4, file.path(figs_dir_path, paste0(cluster_res04,"_padj0.01_DEG.csv")), row.names=FALSE)

# p <- Heatmapplot(seurat.obj, cluster_markers3, cluster_res03)
# ggsave(file.path(figs_dir_path, paste0( cluster_res03,"_heatmap.tiff")),
#       plot = p, units="in", width=size*1.5, height=size*3, dpi=300, compression = 'lzw')

# # make heatmap
# p <- Heatmapplot(seurat.obj, cluster_markers4, cluster_res04)
# ggsave(file.path(figs_dir_path, paste0( cluster_res04,"_heatmap.tiff")),
#       plot = p, units="in", width=size*1.5, height=size*3, dpi=300, compression = 'lzw')

sessionInfo()