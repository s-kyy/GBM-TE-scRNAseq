
# import 2 files, seurat object + cellular marker CSV (one at a time)

library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("At least 3 filepaths must be supplied: [xxx.rds] [xxx_markers_0.3.csv] [xxx_markers_0.4.csv]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1]) & file.exists(args[2]) & file.exists(args[3])) { 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    # parent_dir_path_obj <- dirname(obj_path)
    # parent_dir_name_obj <- basename(parent_dir_path_obj)
    marker_path3 <- args[2]
    marker_path4 <-args[3]
#     parent_dir_path_marker <- dirname(marker_path)
#     parent_dir_name_marker <- basename(parent_dir_path_marker)
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
# library(reshape2)
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
DefaultAssay(seurat.obj) <- "RNA"
# cluster_markers3 <- read.csv(marker_path3)
# cluster_markers4 <- read.csv(marker_path4)
# cluster_res03 <- "integrated_snn_res.0.3"
# cluster_res04 <- "integrated_snn_res.0.4"
size    <- 5

figs_dir_path <- file.path(parent_dir_path_obj, "figs_celltypes")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

print(seurat.obj)

# List of known Markers for each brain cell type
known_markers <- list(
  mes = c('VIM','SOX2','LUM','AXL1'),
  ne = c('SOX2','VIM','NES'),
  rg.pan = c('VIM','NES','SOX2','HES1','EMX1','FGF10','RPS6','NOTCH1','HES5','ATP1A2'),
  rg.early = c('VIM','NES','SOX2','LEF1','GLI2','GFAP'),
  oRG = c('VIM','NES','SOX3','HOPX','FAM107A','IL6ST','PTPRZ1','TNC','LIFR','SOX2','PAX6','SFRP1','ITGB5'),
  vRG = c('VIM','NES','SOX4','FZD8','DAB1'),
  tRG = c('VIM','NES','SOX5','FZD8','TFAP2C','ZIC2','MOXD1','C1ORF61','CRYAB','NFATC2','GPX3'),
  opc = c('OLIG2','OMG','SOX10','SOX2','PDGFRA','DCX'),
  apc = c('AQP4','COL2A1'),
  ol = c('SOX10','MBP','PLP','MAG'),
  ac = c('AQP4','GFAP'),
  npc = c('PAX6','EOMES','SOX2','SFRP1','NEUROG1','NEUROD4','PENK','SSTR2','OXTR','OTX2','FOXG1','LHX2','SIX3','NKX2-1'),
  ipc = c('EOMES','PPP1R17','NKX2-1','ASCL1','DLX1','DLX5','LHX8','LHX6','DLL1','DLL3','HES6','ASCL1','CCND2','NEUROG1','NRN1','STMN2','NEUROD6','DCX','PAX6','SFRP1'),
  ipc.early = c('NHLH1','SLC1A3'),
  ipc.late = c('NHLH2','TP53I11','PPP1R17'),
  n = c('DCX','NEUROD2','GRIA2','GRIA4','SST','RBFOX3','STMN2'),
  n.immature = c('DCX','CXCR4','CXCR7','ERBB4','NRN1','NDST4'),
  n.ex = c('DCX','NEUROD6'),
  n.in = c('DCX','GAD1'),
  inter = c('GAD1', 'PDE4DIP'),
  inter.early = c('CALB2','SST','TAC3','LHX6','DLX6-AS1'),
  inter.late = c('CALB1','CCK','VIP'),
  mg = c('C1QC','TMEM119','P2RY12','CX3CR1','HEXB','TGFB1','CD163'),
  endo = c('PECAM1', 'VWF', 'CLDN5'),
  rbc = c('HBB'),
  macro = c('CX3CR1','PTPRC','SIGLEC1'),
  dend = c('SIGLEC1')
)


#### ======================================================================= ###
#### Create UMAP plots with Known Markers ####
#### ======================================================================= ###

makeKnownMarkerPlots <- function(celltype, markers, obj) {

  DefaultAssay(obj) <- "RNA"

  for (gene in 1:length(markers)) {
    p <- FeaturePlot(obj, features = markers[gene])
    p <- p + theme(axis.line=element_blank(),
            axis.text.x=element_blank(),axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    ggsave(file.path(figs_dir_path, paste0(filename, "_", celltype, "_", markers[gene], "_FeaturePlot.tiff")),
        plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
  }
}

# Filter out genes not expressed in assay
DefaultAssay(seurat.obj) <- "RNA"
condition <- lapply(known_markers, function(x)  x %in% rownames(seurat.obj))
# Create UMAPs per gene
for (i in names(known_markers)) {
  present_markers <- known_markers[[i]][condition[[i]]]
  if (length(present_markers) == 0) {
    print(paste(i, "does not have any matching genes in the assay"))
    next
  } else {
    print(paste(length(present_markers), "of", length(known_markers[[i]]),"genes from celltype",i,"in current assay."))
    print(present_markers)
    makeKnownMarkerPlots(i,present_markers,seurat.obj)
  }
}

#### End of Script #### 
sessionInfo()