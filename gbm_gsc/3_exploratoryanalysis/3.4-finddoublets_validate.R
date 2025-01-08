
# import 2 files, seurat object + cellular marker CSV (one at a time)

# .libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
# .libPaths()
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least 3 filepaths must be supplied: [xxx.rds] [xxx_markers_0.3.csv] [xxx_markers_0.4.csv]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2]) && file.exists(args[3])) { 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    marker_path3 <- args[2]
    marker_path4 <-args[3]
  } else {
    stop("Filepaths provided do not exist. Exiting...", call.=FALSE)
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

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(108)

#### ===================================================================== ####
#### Load Datasets ####
#### ===================================================================== ####
seurat_obj <- readRDS(obj_path) 
DefaultAssay(seurat_obj) <- "RNA"
# cluster_markers3 <- read.csv(marker_path3)
# cluster_markers4 <- read.csv(marker_path4)
# cluster_res03 <- "integrated_snn_res.0.3"
# cluster_res04 <- "integrated_snn_res.0.4"
size    <- 5

figs_dir_path <- file.path(parent_dir_path_obj, "figs_celltypes")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

print(seurat_obj)

# List of known markers for developing brain cells 
stem_markers <- list(
  stem = c('VIM', 'SOX2'),
  mes = c('LUM','AXL1'),
  ne = c('NES'),
  rg.pan = c('HES1','EMX1','FGF10','RPS6','NOTCH1','HES5','ATP1A2'),
  oRG = c('HOPX','FAM107A','IL6ST','PTPRZ1','TNC','LIFR','PAX6','SFRP1','ITGB5'),
  vRG = c('FZD8','DAB1'),
  tRG = c('FZD8','TFAP2C','ZIC2','MOXD1','C1ORF61','CRYAB','NFATC2','GPX3'),
  npc = c('PAX6','EOMES','SOX2','SFRP1','NEUROG1','NEUROD4','PENK','SSTR2','OXTR','OTX2','FOXG1','LHX2','SIX3','NKX2-1'),
  ipc = c('EOMES','PPP1R17','NKX2-1','ASCL1','DLX1','DLX5','LHX8','LHX6','DLL1','DLL3','HES6','ASCL1','CCND2','NEUROG1','NRN1','STMN2','NEUROD6','DCX','PAX6','SFRP1'),
  ipc.early = c('NHLH1','SLC1A3'),
  ipc.late = c('NHLH2','TP53I11','PPP1R17')
)

# List of known Markers for each brain cell type cross-checked with Allen Brain Institute, PangeoDB (developing brain)
known_markers_unique <- list(
  stem = c('VIM', 'SOX2'),
  mes = c('LUM','AXL1'),
  ne = c('NES'),
  opc = c('OMG','SOX10','PDGFRA'),
  ol = c('MBP','PLP','MAG','MOG','OLIG2'),
  ac = c('AQP4', 'GFAP', 'SLC1A2'), # GFAP also expressed by endo, and dev. astroycte (accurate in pangaodb)
  prog = c('PAX6','EOMES','SFRP1','NEUROG1','NKX2-1'), #neural / intermediate progenitor cells
  n.immature = c('NEUROG2'),
  n.early = c('DCX'),
  n.ex = c('SLC17A7'),
  n.in = c('GAD1'),
  inter = c('PVALB'), # PangaoDB denies GAD1 as a marker, no info on PDE4DIP
  mg = c('C1QC','CX3CR1','FYB'),
  pericyte = c('MUSTN1'),
  endo = c('CLDN5','PECAM1', 'VWF'),
  rbc = c('HBB'),
  macro = c('CX3CR1','PTPRC','SIGLEC1'),
  tcell = c('CD3D', 'CD3G'), # 'CD25', 'CD4', 'CD8', 'ZAP70', 'SRK', 'CTLA4', 'FOXP3', 'GITR', 'IKZF2')
  dend = c('SIGLEC1')
)

features_unique <- unlist(known_markers_unique)
cell_types <- names(known_markers_unique)
features_stem <- unlist(stem_markers)
stem_types <- names(stem_markers)

if (grepl("healthy", parent_dir_path_obj, fixed = TRUE)) {
  sample_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#CC79A7", "#ff716e",
    "#999999", "#0072B2", "#194c76", 
    "#D55E00", "#3a4f41", "#6699cc", "#713e5a")
  female_samples <- c("SRR9262922", "SRR9262937",
                        "SRR9264382", "SRR9264383",
                        "SRR9264388")

} else if (grepl("gbm", parent_dir_path_obj, fixed = TRUE)) {
  sample_palette <- c(
    "#E69F00", "#56B4E9", "#009E73", 
    "#F0E442", "#CC79A7", "#ff716e",
    "#999999", "#0072B2", "#194c76")
  female_samples <- "SF11209"
}

#### ======================================================================= ###
#### Create UMAP plots with Known Markers ####
#### ======================================================================= ###

makeUMAPPlot <- function(obj, features) {
  p <- FeaturePlot(obj, features = features)
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

makeKnownMarkerPlots <- function(celltype, markers, obj) {

  DefaultAssay(obj) <- "RNA"

  for (gene in 1:length(markers)) {
    p <- makeUMAPPlot(obj, features = markers[gene])
    ggsave(file.path(figs_dir_path, paste0(celltype, "_", markers[gene], "_FeaturePlot.tiff")),
        plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
  }
}

# Filter out genes not expressed in assay
DefaultAssay(seurat_obj) <- "RNA"

# Create UMAPs per unique marker 
condition <- lapply(known_markers_unique, function(x)  x %in% rownames(seurat_obj))
for (i in names(known_markers_unique)) {
  present_markers <- known_markers_unique[[i]][condition[[i]]]
  if (length(present_markers) == 0) {
    print(paste(i, "does not have any matching genes in the assay"))
    next
  } else {
    print(paste(length(present_markers), "of", length(known_markers_unique[[i]]),"genes from celltype",i,"in current assay."))
    print(present_markers)
    makeKnownMarkerPlots(i,present_markers,seurat_obj)
  }
}

# Create UMAPs per stem cell marker
condition <- lapply(stem_markers, function(x)  x %in% rownames(seurat_obj))
for (i in names(stem_markers)) {
  present_markers <- stem_markers[[i]][condition[[i]]]
  if (length(present_markers) == 0) {
    print(paste(i, "does not have any matching stem cell markers in the assay"))
    next
  } else {
    print(paste(length(present_markers), "of", length(stem_markers[[i]]),"stem cell markers from celltype",i,"in current assay."))
    print(present_markers)
    makeKnownMarkerPlots(i,present_markers,seurat_obj)
  }
}

# check for higher MALAT1 expression in nuclei (healthy). 
# should be expressed evenly everywhere (every cluster, sample etc). 
# https://kb.10xgenomics.com/hc/en-us/articles/360004729092-Why-do-I-see-high-levels-of-Malat1-in-my-gene-expression-data
# https://www.biorxiv.org/content/10.1101/2024.07.14.603469v2

if (("MALAT1") %in% rownames(seurat_obj)) {
  DefaultAssay(seurat_obj) <- "RNA"
  p <- makeUMAPPlot(seurat_obj, features = "MALAT1")
  ggsave(file.path(figs_dir_path, paste0("MALAT1_FeaturePlot.tiff")),
    plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
} else {
  print("MALAT1 does not exist in this seurat_obj RNA assay")
}

#### ======================================================================= ###
#### Create Dot plots with Known Markers ####
#### ======================================================================= ###

# Using scale.data slot for Mean Average Expression

# seurat_obj$`integrated_snn_res.0.3` <- factor(seurat_obj$`integrated_snn_res.0.3`, 
#                                               levels = sort(unique(seurat_obj$`integrated_snn_res.0.3`)))
# seurat_obj$`integrated_snn_res.0.4` <- factor(seurat_obj$`integrated_snn_res.0.4`, 
#                                               levels = sort(unique(seurat_obj$`integrated_snn_res.0.4`)))

DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "integrated_snn_res.0.3"

condition <- features_unique %in% rownames(seurat_obj)
print(paste("Length of unique features to create DotPlot", length(features_unique)))

p <- DotPlot(seurat_obj, features = rev(unique(features_unique[condition]))) + scale_colour_gradient2() + 
      xlab("Known Unique Brain Cell Markers") + ylab("Clusters") + coord_flip() #+ 
      # theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
ggsave(file.path(figs_dir_path, paste0("uniqueMarkers_DotPlot03.tiff")),
    plot = p, units="in", width=size*1.1, height=size*1.3, dpi=300, compression = 'lzw')

Idents(seurat_obj) <- "integrated_snn_res.0.4"

p <- DotPlot(seurat_obj, features = rev(unique(features_unique[condition]))) +   scale_colour_gradient2() + 
      xlab("Known Unique Brain Cell Markers") + ylab("Clusters") + coord_flip() #+ 
      # theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
ggsave(file.path(figs_dir_path, paste0("uniqueMarkers_DotPlot04.tiff")),
    plot = p, units="in", width=size*1.1, height=size*1.3, dpi=300, compression = 'lzw')

#### End of Script #### 
sessionInfo()