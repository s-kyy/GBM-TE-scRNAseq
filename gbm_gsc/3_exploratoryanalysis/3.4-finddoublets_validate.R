
# import 2 files, seurat object + cellular marker CSV (one at a time)

# .libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
# .libPaths()
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("At least 1 filepath must be supplied: [xxx_ge.rds] [figs_dir_name] [cluster_res]", call.=FALSE)
} else {
  # verify filepaths
  if (file.exists(args[1])) { 
    obj_path <- args[1] 
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
  } else {
    stop("Filepaths provided do not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>1) {
    figs_dir_name <- args[2]
  } else {
    figs_dir_name <- "figs_validatecelltypes"
  }

  if (length(args)==3) {
    cluster_res <- args[3]
  } else {
    cluster_res <- "integrated_snn_res.0.4"
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
library(viridis)
library(rlang)

library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(108)

#### ===================================================================== ####
#### Load Datasets ####
#### ===================================================================== ####
seurat_obj <- readRDS(obj_path) 
DefaultAssay(seurat_obj) <- "RNA"
size    <- 5

figs_dir_path <- file.path(parent_dir_path_obj, figs_dir_name)

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

print(seurat_obj)

# List of known markers for developing brain cells 
stem_markers <- list(
  mes = c('LUM','AXL1'),
  npc = c('VIM','NES','SOX2','HES1','NOTCH2'),
  rg.pan = c('HES1','EMX1','FGF10','RPS6','NOTCH1','HES5','ATP1A2'),
  oRG = c('HOPX','TNC','GFAP','EGFR','FAM107A','IL6ST','PTPRZ1','LIFR','PAX6','SFRP1','ITGB5'),
  vRG = c('FZD8','DAB1'),
  tRG = c('FZD8','TFAP2C','ZIC2','MOXD1','C1ORF61','CRYAB','NFATC2','GPX3'),
  npc = c('PAX6','EOMES','SOX2','SFRP1','NEUROG1','NEUROD4','PENK','SSTR2','OXTR','OTX2','FOXG1','LHX2','SIX3','NKX2-1'),
  ipc = c('EOMES','SATB2','BCL11B','PPP1R17','NKX2-1','ASCL1','DLX1','DLX5','LHX8','LHX6','DLL1','DLL3','HES6','CCND2','NEUROG1','NRN1','STMN2','DCX','PAX6','SFRP1'),
  ipc.early = c('NHLH1','SLC1A3'),
  ipc.late = c('NHLH2','TP53I11','PPP1R17')
)

# List of known Markers for each brain cell type cross-checked with Allen Brain Institute, PangeoDB (developing brain)
known_markers_unique <- list(
  mes = c('LUM','AXL1'),
  npc = c('VIM','NES','SOX2','HES1','NOTCH2'),
  opc = c('OMG','SOX10','PDGFRA','EGFR','OLIG2','APOD'),
  ol = c('MBP','PLP','MAG','MOG'),
  ac = c('AQP4', 'GFAP', 'SLC1A2'), # GFAP also expressed by endo, and dev. astroycte (accurate in pangaodb)
  n.ex = c('SLC17A7','NEUROD6'),
  n.in = c('GAD1','NEUROD6'),
  inter = c('PVALB'), # PangaoDB denies GAD1 as a marker, no info on PDE4DIP
  mg = c('C1QC','CX3CR1','FYB'),
  pericyte = c('MUSTN1'),
  endo = c('CLDN5','PECAM1', 'VWF'),
  macro = c('CX3CR1','PTPRC','SIGLEC1'),
  tcell = c('CD3D', 'CD3G'), # 'CD25', 'CD4', 'CD8', 'ZAP70', 'SRK', 'CTLA4', 'FOXP3', 'GITR', 'IKZF2')
  tcell = c('CD27', 'CD79A'), 
  dend = c('SIGLEC1'),
  proliferation = c('MCM6','MCM7','PCNA','MKI67','TOP2A') # G1, S, G2, M
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

# #### ======================================================================= ###
# #### Create UMAP plots with Known Markers ####
# #### ======================================================================= ###

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
        plot = p , units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
  }
}

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
    plot = p , units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
} else {
  print("MALAT1 does not exist in this seurat_obj RNA assay")
}

#### ======================================================================= ###
#### Create Dot plots with Known Markers ####
#### ======================================================================= ###
# Using scale.data slot for Mean Average Expression

DefaultAssay(seurat_obj) <- "RNA"
# Idents(seurat_obj) <- "integrated_snn_res.0.3"

# condition <- features_unique %in% rownames(seurat_obj)
# print(paste("Length of unique features to create DotPlot", length(features_unique)))

# p <- DotPlot(seurat_obj, features = rev(unique(features_unique[condition]))) + scale_colour_gradient2() + 
#       xlab("Known Unique Brain Cell Markers") + ylab("Clusters") + coord_flip() #+ 
#       # theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
# ggsave(file.path(figs_dir_path, paste0("uniqueMarkers_DotPlot03.tiff")),
#     plot = p , units="in", width=size*1.4, height=size*1.3, dpi=300, compression = 'lzw')

Idents(seurat_obj) <- cluster_res
condition <- features_unique %in% rownames(seurat_obj)
print(paste("Length of unique features to create DotPlot", length(features_unique)))

p <- DotPlot(seurat_obj, features = rev(unique(features_unique[condition]))) +   scale_colour_gradient2() + 
      xlab("Known Unique Brain Cell Markers") + ylab("Clusters") + coord_flip() #+ 
      # theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
ggsave(file.path(figs_dir_path, paste0("uniqueMarkers_DotPlot", cluster_res, ".tiff")),
    plot = p , units="in", width=size*1.4, height=size*1.3, dpi=300, compression = 'lzw')

Idents(seurat_obj) <- cluster_res
condition <- features_stem %in% rownames(seurat_obj)
print(paste("Length of unique stem cell markers to create DotPlot", length(unique(features_stem))))

p <- DotPlot(seurat_obj, features = rev(unique(features_stem[condition]))) +   scale_colour_gradient2() + 
      xlab("Known Unique Brain Cell Markers") + ylab("Clusters") + coord_flip() #+ 
      # theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
ggsave(file.path(figs_dir_path, paste0("uniqueStemMarkers_DotPlot", cluster_res, ".tiff")),
    plot = p , units="in", width=size*1.4, height=size*1.5, dpi=300, compression = 'lzw')  

#### ======================================================================= ###
#### Export cells per sample meta.data ####
#### ======================================================================= ###

if (class(seurat_obj@meta.data$mitoRatio) == "list") {
  # seurat_obj@meta.data$mitoRatio <- unlist(seurat_obj@meta.data$mitoRatio)
  meta <- seurat_obj@meta.data
  meta$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")[,1]
  meta$mitoRatio <- meta$mitoRatio / 100
  seurat_obj@meta.data <- meta
} 

# Table of Cell Counts per Sample
df <- seurat_obj@meta.data %>% 
  group_by(sample) %>% 
  summarise( 
    ncells = n(), 
    noveltyscore.avg = mean(log10GenesPerUMI),
    noveltyscore.sd = sd(log10GenesPerUMI),
    ngene.mad = mad(nGene, constant = 1),
    numi.mad = mad(nUMI, constant = 1),
    mitoRatio.mad = mad(mitoRatio, constant = 1),
    mitoRatio.prefilt.mad = mad(mitoRatio_prefilt, constant = 1)
  )

# Count number of detected genes in this assay and avg novelty score
# first determine number of genes (rows) not expressed in any cell (column)
obj.list <- SplitObject(seurat_obj, split.by="sample")
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
write.csv(df, file= file.path(figs_dir_path, paste0(filename,"_samplecounts.csv")), row.names = F)
rm("df")
gc()

#### ======================================================================= ###
#### Export cells per cluster meta.data ####
#### ======================================================================= ###

cluster_summary <- seurat_obj@meta.data %>%
    group_by(!!sym(cluster_res)) %>%
    summarize(count = n(), .groups = 'drop')
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path, paste0("cluster_summary", cluster_res, ".csv")))

cluster_summary <- seurat_obj@meta.data %>%
    group_by(sample, !!sym(cluster_res)) %>%
    summarize(count = n(), .groups = 'drop') %>%
    as.data.frame() %>% ungroup() %>%
    group_by(sample) %>%
    mutate(perc =  round(count/sum(count)*100, 2)) %>% 
    ungroup
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path, paste0("perc_clustersInsamples_summary", cluster_res, ".csv")))

cluster_summary <- seurat_obj@meta.data %>%
    group_by(sample, !!sym(cluster_res)) %>%
    summarize(count = n(), .groups = 'drop') %>%
    as.data.frame() %>% ungroup() %>%
    group_by(!!sym(cluster_res)) %>%
    mutate(perc =  round(count/sum(count)*100, 2)) %>% 
    ungroup
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path, paste0("perc_samplesInClusters_summary", cluster_res, ".csv")))


p <- cluster_summary %>%
    ggplot(aes(x=.data[[cluster_res]], y=perc, fill=sample)) + 
    labs(x = "", y = "Cells (%)", fill="Sample") +
    geom_bar(stat="identity") +
    theme_classic() + 
    scale_fill_manual(values = sample_palette)
    # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave(file.path(figs_dir_path, paste0("barplot_clusterbysample_", cluster_res, ".tiff")),
    plot = p , units="in", width=size*1.1, height=size*0.7, dpi=300, compression = 'lzw')  

#### End of Script #### 
sessionInfo()