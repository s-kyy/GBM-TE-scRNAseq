#!usr/bin/env Rscript

# .libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
# .libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least 2 filepaths must be supplied: [xxx_ge.rds] [yyy_gte.rds] [figs_dir_name] [cluster_res]", call.=FALSE)
} else if (length(args)>=2) {
  
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2])) { 
    obj_path <- args[1] 
    obj_path_2 <- args[2] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_path_obj_2 <- dirname(obj_path_2)
    # parent_dir_name_obj <- basename(parent_dir_path_obj)
  } else {
    stop("Filepaths provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>2) {
    figs_dir_name <- args[3]
  } else {
    figs_dir_name <- "figs_celltypes_annotated"
  }

  if (length(args)==4) {
    # Name of column from integrated dataset 
    cluster_res_num <- as.numeric(args[4])
    cluster_res <- paste0("integrated_snn_res.", cluster_res_num)
    # Custom column name for celltypes in metadata
    cluster_res_symbol <- unlist(strsplit(as.character(cluster_res_num), ".", fixed=TRUE))
    cluster_col_gbm_neftel <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_gbmneftel")
    cluster_col_gbm_tcga <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_gbmtcga")
  } else {
    cluster_res <- "integrated_snn_res.0.4"
    cluster_col_gbm_neftel <- "int04_gbmneftel"
    cluster_col_gbm_tcga <- "int04_gbmtcga"
  }

} else {
    stop("Error in calling R or filepaths. Exiting...", call.=FALSE)
}

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 
library(dplyr)
library(viridis)
library(rlang)
library(reshape2)
library(car)

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(108)
options(warn=1) #print warning messages as they occur

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
seurat_obj_ge <- readRDS(obj_path) 
seurat_obj_gte <- readRDS(obj_path_2) 
size    <- 5
ndims <- 25

# subdir <- file.path(getwd(), paste0(format(Sys.Date(), "%Y%m%d"), "_", sample_name,"_", filename))
subdir <- parent_dir_path_obj
subdir_2 <- parent_dir_path_obj_2

ifelse(!dir.exists(file.path(subdir)),
        dir.create(file.path(subdir),recursive=T),
        "Directory Exists")
ifelse(!dir.exists(file.path(subdir_2)),
        dir.create(file.path(subdir_2),recursive=T),
        "Directory Exists")

figs_dir_path <- file.path(subdir, paste0(figs_dir_name, "_ge"))
figs_dir_path_2 <- file.path(subdir_2, paste0(figs_dir_name, "_gte"))

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")
ifelse(!dir.exists(figs_dir_path_2),
        dir.create(figs_dir_path_2,recursive=T),
        "Directory Exists")

if (grepl("healthy", subdir, fixed = TRUE)) {
  sample_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#CC79A7", "#ff716e",
    "#999999", "#0072B2", "#194c76", 
    "#D55E00", "#3a4f41", "#6699cc", "#713e5a")
  female_samples <- c("SRR9262922", "SRR9262937",
                        "SRR9264382", "SRR9264383",
                        "SRR9264388")

} else if (grepl("gbm", subdir, fixed = TRUE)) {
  sample_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#CC79A7", "#ff716e",
    "#999999", "#0072B2", "#194c76")
  female_samples <- "SF11209"
}

level_neftel <- c(
  "NEFTEL-AC",
  "NEFTEL-OPC",
  "NEFTEL-NPC1",
  # "NEFTEL-NPC2",
  "NEFTEL-NPC1-2",
  "NEFTEL-NPC1-G1.S",
  "NEFTEL-MES1",
  "NEFTEL-MES2",
  "NEFTEL-MES1-2",
  "NEFTEL-G1.S-G2.M",
  "NA"
)

gbmPalette <- c("#FF6B6B", "#56B4E9", "#80DB84",  # CL, MES, PN,
                   "#985F99", "#F4EC7B",  "#1CB086", "#A5B6AD") # CL-MES, CL-PN, MES-PN, NA

colours_neftel <- c(
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

level_tcga <- c(
  "TCGA-PN",
  "TCGA-CL",
  "TCGA-MES",
  "NA"
)

colours_tcga <- c(
  "#80DB84", # "TCGA-PN",
  "#FF6B6B", # "TCGA-CL",
  "#56B4E9", # "TCGA-MES",
  "#A5B6AD" # "NA"
)

#### =========================================== ####
#### Label cells by Neftel GBM subtypes ####
#### =========================================== ####

cluster_neftel <- c(
  '0'  =  "NEFTEL-MES1-2", # "Mesenchymal", 
  '1'  =  "NEFTEL-NPC1", # "OPC Early", 
  '2'  =  "NEFTEL-AC", # "Astrocyte", 
  '3'  =  "NEFTEL-AC", # "Radial Glia 1",
  '4'  =  "NEFTEL-MES2", # "Hypoxic", 
  '5'  =  "NEFTEL-MES1", # "Microglia", 
  '6'  =  "NEFTEL-G1.S-G2.M", # "Cycling 1", 
  '7'  =  "NEFTEL-G1.S-G2.M", # "Mixed Cycling", 
  '8'  =  "NEFTEL-NPC1-2", # "Migrating Neuron", 
  '9'  =  "NEFTEL-NPC1", # "Pre-OPC", 
  '10' =  "NEFTEL-NPC1-G1.S", # "OPC Early Cycling", 
  '11' =  "NA", # "Cycling 2", 
  '12' =  "NA", # "Oligodendrocyte", 
  '13' =  "NA", # "Immune", 
  '14' =  "NA", # "Unknown", 
  '15' =  "NEFTEL-OPC", # "OPC Late", 
  '16' =  "NA", # "Radial Glia 2", 
  '17' =  "NA", # "Endothelia", 
  '18' =  "NA" # "Tumour Endothelia" 
)

meta <- seurat_obj_ge@meta.data
meta[[cluster_col_gbm_neftel]] <- ""
if (is.factor(meta[[cluster_res]])) { meta[[cluster_res]] <- as.character(meta[[cluster_res]]) }

for (i in 1:length(cluster_neftel)) {
  meta[[cluster_col_gbm_neftel]][which(meta[[cluster_res]] == names(cluster_neftel[i]))]  <- cluster_neftel[i]
}

meta[[cluster_col_gbm_neftel]] <- factor(meta[[cluster_col_gbm_neftel]], levels = level_neftel)
meta[[cluster_col_gbm_neftel]] <- droplevels(meta[[cluster_col_gbm_neftel]])

seurat_obj_ge <- AddMetaData(seurat_obj_ge, meta[cluster_col_gbm_neftel], col.name = cluster_col_gbm_neftel)
seurat_obj_gte <- AddMetaData(seurat_obj_gte, meta[cluster_col_gbm_neftel], col.name = cluster_col_gbm_neftel)

#### =========================================== ####
#### Label cells by TCGA GBM subtypes ####
#### =========================================== ####
cluster_tcga <- c(
  '0'  =  "TCGA-MES", # "Mesenchymal", 
  '1'  =  "TCGA-PN", # "OPC Early", 
  '2'  =  "TCGA-CL", # "Astrocyte", 
  '3'  =  "TCGA-CL", # "Radial Glia 1",
  '4'  =  "TCGA-MES", # "Hypoxic", 
  '5'  =  "TCGA-MES", # "Microglia", 
  '6'  =  "TCGA-PN", # "Cycling 1", 
  '7'  =  "TCGA-PN", # "Mixed Cycling", 
  '8'  =  "TCGA-PN", # "Migrating Neuron", 
  '9'  =  "TCGA-PN", # "Pre-OPC", 
  '10' =  "TCGA-PN", # "OPC Early Cycling", 
  '11' =  "NA", # "Cycling 2", 
  '12' =  "NA", # "Oligodendrocyte", 
  '13' =  "NA", # "Immune", 
  '14' =  "NA", # "Unknown", 
  '15' =  "TCGA-PN", # "OPC Late", 
  '16' =  "NA", # "Radial Glia 2", 
  '17' =  "NA", # "Endothelia", 
  '18' =  "NA" # "Tumour Endothelia" 
)

meta[[cluster_col_gbm_tcga]] <- ""
if (is.factor(meta[[cluster_res]])) { meta[[cluster_res]] <- as.character(meta[[cluster_res]]) }

for (i in 1:length(cluster_tcga)) {
  meta[[cluster_col_gbm_tcga]][which(meta[[cluster_res]] == names(cluster_tcga[i]))]  <- cluster_tcga[i]
}

meta[[cluster_col_gbm_tcga]] <- factor(meta[[cluster_col_gbm_tcga]], levels = level_tcga)
meta[[cluster_col_gbm_tcga]] <- droplevels(meta[[cluster_col_gbm_tcga]])

seurat_obj_ge <- AddMetaData(seurat_obj_ge, meta[cluster_col_gbm_tcga], col.name = cluster_col_gbm_tcga)
seurat_obj_gte <- AddMetaData(seurat_obj_gte, meta[cluster_col_gbm_tcga], col.name = cluster_col_gbm_tcga)

#### =========================================== ####
#### Save objects ####
#### =========================================== ####
glimpse(meta)

write.csv(meta, file.path(figs_dir_path, paste0("metadata_ge_cnv_gbm.csv")))
write.csv(meta, file.path(figs_dir_path_2, paste0("metadata_gte_cnv_gbm.csv")))
rm("meta")

# save rds
saveRDS(seurat_obj_ge, file = file.path(subdir, paste0("gbm_ge_celltypes_cnv_gbm.rds")))
saveRDS(seurat_obj_gte, file = file.path(subdir_2, paste0("gbm_gte_celltypes_cnv_gbm.rds")))

# print("Objects saved. Making figures")

#### =========================================== ####
#### UMAPs by Neftel & TCGA GBM subtypes ####
#### =========================================== ####

makeUMAPPlot <- function(obj, ident, labels=FALSE, mycols) {
  DefaultAssay(obj) <- "integrated"
  Idents(obj) <- ident
  p <- DimPlot(obj, reduction = "umap", group.by = ident, label=labels, label.size = 3.5, cols=mycols) 
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

p <- makeUMAPPlot(obj = seurat_obj_ge, ident = cluster_col_gbm_neftel, mycols=colours_neftel)
ggsave(file.path(figs_dir_path, paste0("cluster_DimPlot_defaultcols_nolabels",cluster_col_gbm_neftel,".tiff")),
  plot = p, units="in", width=size*1.3, height=size*1.2, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj = seurat_obj_ge, ident = cluster_col_gbm_neftel, labels = TRUE, mycols=colours_neftel)
ggsave(file.path(figs_dir_path, paste0("cluster_DimPlot_defaultcols",cluster_col_gbm_neftel,".tiff")),
  plot = p, units="in", width=size*1.3, height=size*1.2, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj = seurat_obj_ge, ident = cluster_col_gbm_tcga, mycols=colours_tcga)
ggsave(file.path(figs_dir_path, paste0("cluster_DimPlot_defaultcols_nolabels",cluster_col_gbm_tcga,".tiff")),
  plot = p, units="in", width=size*1.2, height=size*1.2, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj = seurat_obj_ge, ident = cluster_col_gbm_tcga, labels = TRUE, mycols=colours_tcga)
ggsave(file.path(figs_dir_path, paste0("cluster_DimPlot_defaultcols",cluster_col_gbm_tcga,".tiff")),
  plot = p, units="in", width=size*1.2, height=size*1.2, dpi=300, compression = 'lzw')

p <- makeUMAPPlot(obj = seurat_obj_gte, ident = cluster_col_gbm_neftel, mycols=colours_neftel)
ggsave(file.path(figs_dir_path_2, paste0("cluster_DimPlot_defaultcols_nolabels",cluster_col_gbm_neftel,".tiff")),
  plot = p, units="in", width=size*1.3, height=size*1.2, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj = seurat_obj_gte, ident = cluster_col_gbm_neftel, labels = TRUE, mycols=colours_neftel)
ggsave(file.path(figs_dir_path_2, paste0("cluster_DimPlot_defaultcols",cluster_col_gbm_neftel,".tiff")),
  plot = p, units="in", width=size*1.3, height=size*1.2, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj = seurat_obj_gte, ident = cluster_col_gbm_tcga, mycols=colours_tcga)
ggsave(file.path(figs_dir_path_2, paste0("cluster_DimPlot_defaultcols_nolabels",cluster_col_gbm_tcga,".tiff")),
  plot = p, units="in", width=size*1.2, height=size*1.2, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj = seurat_obj_gte, ident = cluster_col_gbm_tcga, labels = TRUE, mycols=colours_tcga)
ggsave(file.path(figs_dir_path_2, paste0("cluster_DimPlot_defaultcols",cluster_col_gbm_tcga,".tiff")),
  plot = p, units="in", width=size*1.2, height=size*1.2, dpi=300, compression = 'lzw')

#### =========================================== ####
#### Dot Plot by Neftel & TCGA GBM subtypes ####
#### =========================================== ####
# Tips: https://github.com/ycl6/StackedVlnPlot 

DefaultAssay(seurat_obj_ge) <- "RNA"
DefaultAssay(seurat_obj_gte) <- "RNA"

# TCGA
tcga_features <- c('ASCL1', 'OLIG2', 'EGFR', 'GFAP', "TGFB1","RELB") 

condition <- tcga_features %in% rownames(seurat_obj_ge)
print(paste("Length of unique features to create VlnPlot", length(unique(tcga_features[condition]))))
Idents(seurat_obj_ge) <- cluster_col_gbm_tcga
p <- VlnPlot(seurat_obj_ge, features = unique(tcga_features[condition]),
  stack = TRUE, flip=TRUE) + xlab("") + ylab("")  + 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1), legend.position = "none")
ggsave(file.path(figs_dir_path, paste0("vlnplot_markers_",cluster_col_gbm_tcga,".tiff")), 
  plot = p, units="in", width=size*0.5, height=size*0.7, dpi=300, compression = 'lzw')

condition <- tcga_features %in% rownames(seurat_obj_gte)
print(paste("Length of unique features to create VlnPlot", length(unique(tcga_features[condition]))))
Idents(seurat_obj_gte) <- cluster_col_gbm_tcga
p <- VlnPlot(seurat_obj_gte, features = unique(tcga_features[condition]),
  stack = TRUE, flip=TRUE) + xlab("") + ylab("")  + 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1), legend.position = "none")
  # RotatedAxis() + xlab('') + ylab('')
ggsave(file.path(figs_dir_path_2, paste0("vlnplot_markers_",cluster_col_gbm_tcga,".tiff")), 
  plot = p, units="in", width=size*0.5, height=size*0.7, dpi=300, compression = 'lzw')

# NEFTEL
neftel_features <- c(
  "S100B", "GFAP", "SLC1A3", "GLAST", "MLC1", "HOPX", # AC + RG
  "OMG", "PLP1", "PLLP", "TNR", "ALCAM", # OPC (+OLIG1)
  "SOX4", "SOX11", "DCX", # NPC
  "OLIG1", "TNR", # NPC-1
  "STMN1", "STMN2", "STMN4", "DLX5-AS1", "DLX6-AS1", # NPC-2
  "VIM", "ANXA1", "ANXA2", "CHI3L1", "CD44", # MES1
  "ADM", "LDHA", "HILPDA", "ENO2", "DDIT3" # MES2
)

condition <- neftel_features %in% rownames(seurat_obj_ge)
print(paste("Length of unique features to create VlnPlot", length(unique(neftel_features[condition]))))
Idents(seurat_obj_ge) <- cluster_col_gbm_neftel
p <- VlnPlot(seurat_obj_ge, features = unique(neftel_features[condition]),
  stack = TRUE, flip=TRUE) + xlab("") + ylab("") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1), legend.position = "none")
ggsave(file.path(figs_dir_path, paste0("vlnplot_markers_",cluster_col_gbm_neftel,".tiff")), 
  plot = p, units="in", width=size*1.1, height=size*1.5, dpi=300, compression = 'lzw')

condition <- neftel_features %in% rownames(seurat_obj_gte)
print(paste("Length of unique features to create VlnPlot", length(unique(neftel_features[condition]))))
Idents(seurat_obj_gte) <- cluster_col_gbm_neftel
p <- VlnPlot(seurat_obj_gte, features = unique(neftel_features[condition]),
  stack = TRUE, flip=TRUE) + xlab("") + ylab("") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1), legend.position = "none")
ggsave(file.path(figs_dir_path_2, paste0("vlnplot_markers_",cluster_col_gbm_neftel,".tiff")), 
  plot = p, units="in", width=size*1.1, height=size*1.5, dpi=300, compression = 'lzw')

### End of Script
sessionInfo()
