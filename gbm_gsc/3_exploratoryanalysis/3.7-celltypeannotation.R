#!usr/bin/env Rscript

# .libPaths(c("~/scratch/tcga-gbm-R4-lib/x86_64-pc-linux-gnu", .libPaths()))
# .libPaths()
#### Parse Arguments ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least 2 filepaths must be supplied: [xxx_ge.rds] [yyy_gte.rds] [cc_markers.csv] [figs_dir_name] [cluster_res]", call.=FALSE)
} else if (length(args)>=3) {
  
  # verify filepaths
  if (file.exists(args[1]) && file.exists(args[2])) { 
    obj_path <- args[1] 
    obj_path_2 <- args[2] 
    # cc.path <- args[3]
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path_obj <- dirname(obj_path)
    parent_dir_path_obj_2 <- dirname(obj_path_2)
    # parent_dir_name_obj <- basename(parent_dir_path_obj)
    ccgenes_path <- args[3]
  } else {
    stop("Filepaths provided does not exist. Exiting...", call.=FALSE)
  }

  if (length(args)>3) {
    figs_dir_name <- args[4]
  } else {
    figs_dir_name <- "figs_celltypes_annotated"
  }

  if (length(args)==5) {
    # Name of column from integrated dataset 
    cluster_res_num <- as.numeric(args[5])
    cluster_res <- paste0("integrated_snn_res.", cluster_res_num)
    # Custom column name for celltypes in metadata
    cluster_res_symbol <- unlist(strsplit(as.character(cluster_res_num), ".", fixed=TRUE))
    cluster_col <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_celltypes")
    cluster_col_gsc <- paste0("int",cluster_res_symbol[1],cluster_res_symbol[2], "_gsctypes")
    cluster_col_gsc_cc <- paste0(cluster_col_gsc,"_cc")
  } else {
    cluster_res <- "integrated_snn_res.0.4"
    cluster_col <- "int04_celltypes"
    cluster_col_gsc <- "int04_gsctypes"
    cluster_col_gsc_cc <- "int04_gsctypes_cc"
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
cc.genes <- read.csv(ccgenes_path)
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

level_celltypes <- c(
  # "NPC (+TNC)",
  # "NPC (-TNC)",
  "Mesenchymal",
  "Hypoxic",
  "Ex. Neuron",
  "Ex. Neuron-NEUROD6",
  "Ex. Neuron-SATB2",
  "Ex. Neuron-Signaling",
  "In. Neuron",
  "In. Neuron-TP53I11",
  "In. Neuron-BCL11B",
  "PVALB Interneuron",
  "Astrocyte", 
  "Radial Glia 1",
  "Radial Glia 2",
  # "Stressed Astrocyte",
  "Outer Radial Glia",
  "Outer Radial Glia Cycling",
  "Outer Radial Glia Cycling S.G2-M",
  "Outer Radial Glia Cycling S",
  "Outer Radial Glia Cycling G2-M",
  # "Outer Radial Glia S",
  # "Outer Radial Glia G2/M",
  "Pre-OPC", 
  "OPC Early",
  # "OPC S Phase",
  "OPC Early Cycling",
  "OPC Early Cycling S",
  "OPC Late",
  "Oligodendrocyte", 
  "Microglia", 
  "Immune",
  "Migrating Neuron",
  "Endothelia",
  "Tumour Endothelia",
  "Purkinje",
  "Cycling 1",
  "Cycling 1 S.G2-M",
  "Cycling 2",
  "Cycling 2 G2-M",
  "Mixed Cycling",
  "Mixed Cycling S.G2-M",
  # "Cycling",
  # "Cycling G1.S",
  # "Cycling (uncontrolled)",
  # "Cycling G2/M",
  "G2/M Phase",
  "Dying Cell", # high MALAT1
  "Unknown"
)

known_markers <- c()

#### =========================================== ####
#### Add cell type annotations ####
#### =========================================== ####

DefaultAssay(seurat_obj_ge) <- "RNA"
DefaultAssay(seurat_obj_gte) <- "RNA"

# Healthy
if (grepl("healthy", subdir, fixed = TRUE)) {

  known_markers <- c(
    'SLC17A7','STMN2','NRN1',# excitatory neurons 
    'NEUROD6', #Ex. neuron (NeuroD6) 
    'SATB2', # Ex. neuron (NeuroD6) 
    'ENC1','TESPA1','GRIA2','GRIA4','HOPX', # Signaling neurons (+MEG3,+FAM153CP )
    'GAD1','DLX1','LHX6', # inhibitory neurons 
    'TP53I11', # in. neuron (TP53I11)
    'BCL11B', # in. neuron (BCL11B)
    'PVALB', # PV.interneurons / GABA-ergic interneurons + GAD1
    'AQP4', 'GFAP', 'SLC1A2','SLC1A3','ATP1A2', # astrocytes
    'OLIG2','PDGFRA','OMG','SOX10','PTPRZ1', # OPC Late
    'APOD','MBP', 'MAG', 'MOG','CRYAB', # oligodendrocytes
    'C1QC', 'CX3CR1','PTPRC', # microglia
    'CLDN5', 'PECAM1', 'VWF','MUSTN1','VIM','RPS6', # endothelial
    'PCP4' # Purkinje (motor neuron)
  )

  # GE + save rds
  celltypes <- c(
    '0' = "Oligodendrocyte",
    '1' = "Astrocyte", 
    '2' = "Ex. Neuron",
    '3' = "OPC Late",
    '4' = "Ex. Neuron-Signaling",#"Progenitor",
    '5' = "Ex. Neuron-Signaling",#"Progenitor",
    '6' = "Ex. Neuron",
    '7' = "In. Neuron",
    '8' = "Ex. Neuron-NEUROD6",
    '9' = "In. Neuron-TP53I11",
    '10' = "PVALB Interneuron",
    '11' = "Ex. Neuron-NEUROD6",
    '12' = "Ex. Neuron-NEUROD6",
    '13' = "Ex. Neuron",
    '14' = "Endothelia",
    '15' = "Ex. Neuron-NEUROD6",
    '16' = "Ex. Neuron-NEUROD6",
    '17' = "In. Neuron-TP53I11",
    '18' = "In. Neuron-BCL11B",
    '19' = "Ex. Neuron-NEUROD6",
    '20' = "Microglia",
    '21' = "Ex. Neuron-SATB2",
    '22' = "PVALB Interneuron",
    '23' = "Purkinje",
    '24' = "PVALB Interneuron",
    '25' = "Oligodendrocyte"
  )
  print(known_markers)
  print(celltypes)
  meta <- seurat_obj_ge@meta.data
  meta[[cluster_col]] <- ""
  if (is.factor(meta[[cluster_res]])) { meta[[cluster_res]] <- as.character(meta[[cluster_res]]) }

  for (i in 1:length(celltypes)) {
    meta[[cluster_col]][which(meta[[cluster_res]] == names(celltypes[i]))]  <- celltypes[i]
  }

  meta[[cluster_col]] <- factor(meta[[cluster_col]], levels = level_celltypes)
  meta[[cluster_col]] <- droplevels(meta[[cluster_col]])

  seurat_obj_ge <- AddMetaData(seurat_obj_ge, 
                                meta[cluster_col], 
                                col.name = cluster_col)

  # GTE
  seurat_obj_gte <- AddMetaData(seurat_obj_gte, 
                                meta[cluster_col], 
                                col.name = cluster_col)
  # save metadata
  if (class(meta$mitoRatio) == "data.frame") { 
    # meta$mitoRatio <- meta$mitoRatio[,1] 
    meta$mitoRatio <- NULL
  } else if (class(meta$mitoRatio) == "list") {
    # meta$mitoRatio <- meta$mitoRatio[[1]]
    meta$mitoRatio <- NULL
  }
  meta <- meta %>% dplyr::mutate(across(everything(), as.character))
  write.csv(meta, file.path(figs_dir_path, paste0("metadata_ge.csv")))
  
  meta <- seurat_obj_gte@meta.data
  meta <- meta %>% dplyr::mutate(across(everything(), as.character))
  if (class(meta$mitoRatio) == "data.frame") { 
    # meta$mitoRatio <- meta$mitoRatio[,1] 
    meta$mitoRatio <- NULL
  } else if (class(meta$mitoRatio) == "list") {
    # meta$mitoRatio <- meta$mitoRatio[[1]]
    meta$mitoRatio <- NULL
  }
  write.csv(meta, file.path(figs_dir_path_2, paste0("metadata_gte.csv")))
  rm("meta")

  # save rds
  saveRDS(seurat_obj_ge, file = file.path(subdir, paste0("healthy_ge_celltypes.rds")))
  saveRDS(seurat_obj_gte, file = file.path(subdir_2, paste0("healthy_gte_celltypes.rds")))

# GBM
} else if (grepl("gbm", subdir, fixed = TRUE)) {

  known_markers <- c(
    'VIM','EGFR', # NPC
    'NRN1','HILPDA','VEGFA','NDRG1', # HYPOXIC/ NEURONAL STRESS 
    'HOPX','TNC','FAM107A','LIFR','PTPRZ1','FABP7','IL6ST', # RG
    'AQP4','GFAP','SLC1A2','SLC1A3', # ASTROCYTE
    'OLIG1','OLIG2','SOX2', # PRE-OPC
    'PDGFRA','NES','HES6','DLL3', # OPC-EARLY : OLIG2
    'OMG','SOX10', # OPC-LATE
    'APOD','MBP','MAG','MOG','CRYAB', # OL
    'C1QC', 'CX3CR1','PTPRC', # MG
    'CD27', 'CD79A', 'IGKC', 'IGHM', #B-cell
    'CD3D', 'CD3G', # T-cell
    'STMN2','DCX','PAX6', # 'NEUROD6', # MATURING/DEVELOPING NEURONS
    'CLDN5', 'PECAM1', 'VWF', # ENDOTHELIAL
    'CD248','LUM', 'RGS5','PDGFRB', # ENDO. TUMOUR
    'MCM7', 'PCNA','MKI67','TOP2A', # proliferating cells
    'MALAT1' #Unknown
  )
  
  # GE
  celltypes <- c(
    '0'  = "Mesenchymal", # "NPC (+TNC)",
    '1'  = "OPC Early", # "OPC Early",
    '2'  = "Astrocyte", # "Astrocyte", 
    '3'  = "Radial Glia 1", # "Astrocyte",
    '4'  = "Hypoxic", # "NPC (-TNC)",#"Neurons (hypoxic)",
    '5'  = "Microglia", # "Microglia", 
    '6'  = "Cycling 1", # "Cycling G1.S", 
    '7'  = "Mixed Cycling", # "Cycling (uncontrolled)", 
    '8'  = "Migrating Neuron", # "Neuron", 
    '9'  = "Pre-OPC", # "Pre-OPC", 
    '10' = "OPC Early Cycling", # "OPC G1.S",
    '11' = "Cycling 2", # "Cycling G2/M",
    '12' = "Oligodendrocyte", # "Oligodendrocyte", 
    '13' = "Immune", # "Immune",
    '14' = "Unknown", # "Unknown", # high MALAT1
    '15' = "OPC Late", # "OPC Late",
    '16' = "Radial Glia 2", # "Stressed Astrocyte", # with stress genes, less FAM107A
    '17' = "Endothelia", # "Endothelia",
    '18' = "Tumour Endothelia" # "Tumour Endothelia"
  )
  print(known_markers)
  print(celltypes)
  meta <- seurat_obj_ge@meta.data
  meta[[cluster_col]] <- ""
  if (is.factor(meta[[cluster_res]])) { meta[[cluster_res]] <- as.character(meta[[cluster_res]]) }

  for (i in 1:length(celltypes)) {
    meta[[cluster_col]][which(meta[[cluster_res]] == names(celltypes[i]))]  <- celltypes[i]
  }

  meta[[cluster_col]] <- factor(meta[[cluster_col]], levels = level_celltypes)
  meta[[cluster_col]] <- droplevels(meta[[cluster_col]])

  seurat_obj_ge <- AddMetaData(seurat_obj_ge, meta[cluster_col], col.name = cluster_col)
  seurat_obj_gte <- AddMetaData(seurat_obj_gte, meta[cluster_col], col.name = cluster_col)

  # GSC metadata
  # Seurat subsetting defaults to "data" slot : https://github.com/satijalab/seurat/issues/4658#issuecomment-868716116
  gsc <- subset(x = seurat_obj_ge, subset = 
    # Filter criteria from Bhaduri et al., 2020
    ((PROM1  > 0) | (FUT4   > 0) | (L1CAM  > 0)) & (SOX2   > 0) & (TLR4   <= 0) & # CD133(PROM1)
    # Filter Immune Cells (T NK B cells)
    (PTPRC <= 0) & (CD27 <= 0) & (CD79A <= 0)
  )
  
  org <- subset(x = gsc, subset = 
    # Filter criteria from Bhaduri et al., 2020 (# Also upregulated : (TNC > 0) & (LIFR > 0) & (HOPX > 0) & (FAM107A > 0) & (IL6ST > 0) )
    (PTPRZ1 > 0) & 
    # Filter criteria from Wang et al., 2020
    (SLC1A3  > 0) &  # GLAST
    (FABP7  > 0) &   # BLBP
    (HOPX   > 0) 
  )

  # Label GSC cells in metadata
  meta <- seurat_obj_ge@meta.data
  meta[['gsc']] <- 0
  
  print(paste("Dimensions of gsc", dim(gsc)[2]))
  barcode_gsc <- colnames(gsc)
  meta[barcode_gsc,'gsc'] <- 1
  table(meta[,'gsc']) # preview labels for GSC vs non-GSC cells

  # Label oRG cells in metadata
  meta[['oRG']] <- 0
  
  print(paste("Dimensions of oRG", dim(org)[2]))
  barcode_org <- colnames(org)
  meta[barcode_org,'oRG'] <- 1
  table(meta[,'oRG']) # preview labels for oRG vs non-oRG cells

  seurat_obj_ge <- AddMetaData(seurat_obj_ge, meta['gsc'], col.name = 'gsc')
  seurat_obj_ge <- AddMetaData(seurat_obj_ge, meta['oRG'], col.name = 'oRG')
  seurat_obj_gte <- AddMetaData(seurat_obj_gte, meta['gsc'], col.name = 'gsc')
  seurat_obj_gte <- AddMetaData(seurat_obj_gte, meta['oRG'], col.name = 'oRG')
  
  # Compute Cell cycling scores for both objects
  s.genes <- unique(cc.genes$s.genes[nzchar(cc.genes$s.genes)])       # extract cc related genes from each column
  g2m.genes <- unique(cc.genes$g2m.genes[nzchar(cc.genes$g2m.genes)]) # extract cc related genes from each column

  print(paste("Number of S-phase related genes", length(s.genes)))
  print(s.genes)
  print(paste("Number of G2/M-phase related genes", length(g2m.genes)))
  print(g2m.genes)

  seurat_obj_ge <- CellCycleScoring(seurat_obj_ge, 
                                    s.features = s.genes, 
                                    g2m.features = g2m.genes, 
                                    set.ident = FALSE)
  seurat_obj_gte@meta.data$S.Score <- seurat_obj_ge@meta.data$S.Score
  seurat_obj_gte@meta.data$G2M.Score <- seurat_obj_ge@meta.data$G2M.Score
  seurat_obj_gte@meta.data$Phase <- seurat_obj_ge@meta.data$Phase
  prop.table(table(seurat_obj_ge@meta.data[c(cluster_res,'Phase')]),margin=1)*100

  # Label cluster_col_gsc & cluster_col_gsc_cc
  meta <- seurat_obj_ge@meta.data
  if (is.factor(meta[[cluster_col]])) { meta[[cluster_col]] <- as.character(meta[[cluster_col]]) }
  if (is.factor(meta[[cluster_res]])) { meta[[cluster_res]] <- as.character(meta[[cluster_res]]) }
  meta[[cluster_col_gsc]] <- meta[[cluster_col]]
  meta[[cluster_col_gsc_cc]] <- meta[[cluster_col]]

  barcode_org <- colnames(org)
  for (cell in barcode_org){
    if (meta[cell, "Phase"] == "S" || meta[cell, "Phase"] == "G2M") {
      meta[cell, cluster_col_gsc] <- "Outer Radial Glia Cycling" 
    } else {
      meta[cell, cluster_col_gsc] <- "Outer Radial Glia"
    }
  }

  meta[[cluster_col_gsc]] <- factor(meta[[cluster_col_gsc]], levels = level_celltypes)
  meta[[cluster_col_gsc]] <- droplevels(meta[[cluster_col_gsc]])

  # Label Cycling oRGs & GSCs
  # meta <- seurat_obj_ge@meta.data
  if (is.factor(meta[[cluster_col]])) { meta[[cluster_col]] <- as.character(meta[[cluster_col]]) }
  if (is.factor(meta[[cluster_res]])) { meta[[cluster_res]] <- as.character(meta[[cluster_res]]) }
  meta[[cluster_col_gsc_cc]] <- meta[[cluster_col]]

  meta[[cluster_col_gsc_cc]][which(meta[[cluster_res]] == '6' )]  <-  "Cycling 1 S.G2-M" # "Cycling"
  meta[[cluster_col_gsc_cc]][which(meta[[cluster_res]] == '7' )]  <-  "Mixed Cycling S.G2-M" # "Cycling G2/M"
  meta[[cluster_col_gsc_cc]][which(meta[[cluster_res]] == '10' )]  <- "OPC Early Cycling S" # "OPC S Phase"
  meta[[cluster_col_gsc_cc]][which(meta[[cluster_res]] == '11' )]  <- "Cycling 2 G2-M" # "G2/M Phase"

  barcode_org <- colnames(org)
  for (cell in barcode_org){
    if (meta[cell, "Phase"] == "S" && meta[cell, "Phase"] == "G2M") {
      meta[cell, cluster_col_gsc_cc] <- "Outer Radial Glia Cycling S.G2-M" 
    } else if (meta[cell, "Phase"] == "S") {
      meta[cell, cluster_col_gsc_cc] <- "Outer Radial Glia Cycling S"
    } else if (meta[cell, "Phase"] == "G2M") {
      meta[cell, cluster_col_gsc_cc] <- "Outer Radial Glia Cycling G2-M"
    } else {
      meta[cell, cluster_col_gsc_cc] <- "Outer Radial Glia"
    }
  }
  
  meta[[cluster_col_gsc_cc]] <- factor(meta[[cluster_col_gsc_cc]], levels = level_celltypes)
  meta[[cluster_col_gsc_cc]] <- droplevels(meta[[cluster_col_gsc_cc]])

  seurat_obj_ge <- AddMetaData(seurat_obj_ge, meta[cluster_col_gsc], col.name = cluster_col_gsc)
  seurat_obj_gte <- AddMetaData(seurat_obj_gte, meta[cluster_col_gsc], col.name = cluster_col_gsc)
  seurat_obj_ge <- AddMetaData(seurat_obj_ge, meta[cluster_col_gsc_cc], col.name = cluster_col_gsc_cc)
  seurat_obj_gte <- AddMetaData(seurat_obj_gte, meta[cluster_col_gsc_cc], col.name = cluster_col_gsc_cc)

  glimpse(seurat_obj_ge@meta.data)
  glimpse(seurat_obj_gte@meta.data)

  if (class(meta$mitoRatio) == "data.frame") { 
    # meta$mitoRatio <- meta$mitoRatio[,1] 
    meta$mitoRatio <- NULL
  } else if (class(meta$mitoRatio) == "list") {
    # meta$mitoRatio <- meta$mitoRatio[[1]]
    meta$mitoRatio <- NULL
  }
  meta <- meta %>% dplyr::mutate(across(everything(), as.character))
  write.csv(meta, file.path(figs_dir_path, paste0("metadata_ge.csv")))
  
  meta <- seurat_obj_gte@meta.data
  meta <- meta %>% dplyr::mutate(across(everything(), as.character))
  if (class(meta$mitoRatio) == "data.frame") { 
    # meta$mitoRatio <- meta$mitoRatio[,1] 
    meta$mitoRatio <- NULL
  } else if (class(meta$mitoRatio) == "list") {
    # meta$mitoRatio <- meta$mitoRatio[[1]]
    meta$mitoRatio <- NULL
  }
  write.csv(meta, file.path(figs_dir_path_2, paste0("metadata_gte.csv")))
  rm("meta")

  # save rds
  if (class(seurat_obj_ge@meta.data$mitoRatio) == "data.frame" || class(seurat_obj_ge@meta.data$mitoRatio) == "list") { seurat_obj_ge@meta.data$mitoRatio$mitoRatio <- NULL } 
  if (class(seurat_obj_gte@meta.data$mitoRatio) == "data.frame"|| class(seurat_obj_gte@meta.data$mitoRatio) == "list") { seurat_obj_gte@meta.data$mitoRatio$mitoRatio <- NULL } 
  saveRDS(seurat_obj_ge, file = file.path(subdir, paste0("gbm_ge_celltypes.rds")))
  saveRDS(seurat_obj_gte, file = file.path(subdir_2, paste0("gbm_gte_celltypes.rds")))
}


#### =========================================== ####
#### Dot Plots with known markers ####
#### =========================================== ####

DefaultAssay(seurat_obj_ge) <- "RNA"
Idents(seurat_obj_ge) <- cluster_col

condition <- known_markers %in% rownames(seurat_obj_ge)
print(paste("Length of unique features to create DotPlot", length(unique(known_markers))))

p <- DotPlot(seurat_obj_ge, features = rev(unique(known_markers[condition]))) + scale_colour_gradient2() + 
      xlab("") + ylab("Cell Types") + coord_flip() + 
      theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
ggsave(file.path(figs_dir_path, paste0("celltypeXuniqueMarkers_DotPlot_",cluster_col,".tiff")),
    plot = p, units="in", width=size*1.6, height=size*1.8, dpi=300, compression = 'lzw')

DefaultAssay(seurat_obj_gte) <- "RNA"
Idents(seurat_obj_gte) <- cluster_col

condition <- known_markers %in% rownames(seurat_obj_gte)
print(paste("Length of unique features to create DotPlot", length(unique(known_markers))))

p <- DotPlot(seurat_obj_gte, features = rev(unique(known_markers[condition]))) + scale_colour_gradient2() + 
      xlab("") + ylab("Cell Types") + coord_flip() + 
      theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
ggsave(file.path(figs_dir_path_2, paste0("celltypeXuniqueMarkers_DotPlot_",cluster_col,".tiff")),
    plot = p, units="in", width=size*1.6, height=size*1.8, dpi=300, compression = 'lzw')


#### =========================================== ####
#### cluster celltype summary with known markers ####
#### =========================================== ####

## Export table counting cells per cluster across all samples
cluster_summary <- seurat_obj_ge@meta.data %>%
    group_by(!!sym(cluster_res)) %>%
    summarize(count = n(), .groups = 'drop')
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path, paste0("cluster_summary_",cluster_res,".csv")))

cluster_summary <- seurat_obj_gte@meta.data %>%
    group_by(!!sym(cluster_res)) %>%
    summarize(count = n(), .groups = 'drop')
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path_2, paste0("cluster_summary_",cluster_res,".csv")))

## Export table counting cells per celltype across all samples
cluster_summary <- seurat_obj_ge@meta.data %>%
    group_by(!!sym(cluster_col)) %>%
    summarize(count = n(), .groups = 'drop')
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path, paste0("celltype_summary_",cluster_col,".csv")))

cluster_summary <- seurat_obj_gte@meta.data %>%
    group_by(!!sym(cluster_col)) %>%
    summarize(count = n(), .groups = 'drop')
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path_2, paste0("celltype_summary_",cluster_col,".csv")))

## Export table counting cells per celltype in each sample
cluster_summary <- seurat_obj_ge@meta.data %>%
    group_by(sample, !!sym(cluster_col)) %>%
    summarize(count = n(), .groups = 'drop') %>%
    as.data.frame() %>% ungroup() %>%
    group_by(sample) %>%
    mutate(perc =  round(count/sum(count)*100, 2)) #%>% ungroup()
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path, paste0("perc_clustersInsamples_summary_",cluster_col,".csv")))

cluster_summary <- seurat_obj_gte@meta.data %>%
    group_by(sample, !!sym(cluster_col)) %>%
    summarize(count = n(), .groups = 'drop') %>%
    as.data.frame() %>% ungroup() %>%
    group_by(sample) %>%
    mutate(perc =  round(count/sum(count)*100, 2)) #%>% ungroup()
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path_2, paste0("perc_clustersInsamples_summary_",cluster_col,".csv")))

## Export table counting cells per sample in each celltype + Bar plot with cell counts per cell type grouped by sample
cluster_summary <- seurat_obj_ge@meta.data %>%
    group_by(sample, !!sym(cluster_col)) %>%
    summarize(count = n(), .groups = 'drop') %>%
    as.data.frame() %>% ungroup() %>%
    group_by(!!sym(cluster_col)) %>%
    mutate(perc =  round(count/sum(count)*100, 2)) #%>% ungroup()
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path, paste0("perc_samplesInCelltype_summary_",cluster_col,".csv")))

p <- cluster_summary %>%
    ggplot(aes(x=!!sym(cluster_col), y=perc, fill=sample)) + 
    labs(x = "", y = "Cells (%)", fill="Sample") +
    geom_bar(stat="identity") +
    theme_classic() + 
    scale_fill_manual(values = sample_palette)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave(file.path(figs_dir_path, paste0("barplot_celltypeBysample.tiff")),
    plot = p, units="in", width=size*1.1, height=size*0.7, dpi=300, compression = 'lzw')  


cluster_summary <- seurat_obj_gte@meta.data %>%
    group_by(sample, !!sym(cluster_col)) %>%
    summarize(count = n(), .groups = 'drop') %>%
    as.data.frame() %>% ungroup() %>%
    group_by(!!sym(cluster_col)) %>%
    mutate(perc =  round(count/sum(count)*100, 2)) #%>% ungroup()
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path_2, paste0("perc_samplesInCelltype_summary_",cluster_col,".csv")))

p <- cluster_summary %>%
    ggplot(aes(x=.data[[cluster_col]], y=perc, fill=sample)) + 
    labs(x = "", y = "Cells (%)", fill="Sample") +
    geom_bar(stat="identity") +
    theme_classic() + 
    scale_fill_manual(values = sample_palette)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave(file.path(figs_dir_path_2, paste0("barplot_celltypeBysample_",cluster_col,".tiff")),
    plot = p, units="in", width=size*1.1, height=size*0.7, dpi=300, compression = 'lzw')  

## Export table counting cells per sample in each cluster + Bar plot with cell counts per cell type grouped by cluster
cluster_summary <- seurat_obj_ge@meta.data %>%
    group_by(sample, !!sym(cluster_res)) %>%
    summarize(count = n(), .groups = 'drop') %>%
    as.data.frame() %>% ungroup() %>%
    group_by(!!sym(cluster_res)) %>%
    mutate(perc =  round(count/sum(count)*100, 2)) #%>% ungroup()
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path, paste0("perc_samplesInCluster_summary_",cluster_res,".csv")))

p <- cluster_summary %>%
    ggplot(aes(x=.data[[cluster_res]], y=perc, fill=sample)) + 
    labs(x = "", y = "Cells (%)", fill="Sample") +
    geom_bar(stat="identity") +
    theme_classic() + 
    scale_fill_manual(values = sample_palette)
    # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave(file.path(figs_dir_path, paste0("barplot_clusterbysample_",cluster_res,".tiff")),
    plot = p, units="in", width=size*1.1, height=size*0.7, dpi=300, compression = 'lzw')  


cluster_summary <- seurat_obj_gte@meta.data %>%
    group_by(sample, !!sym(cluster_res)) %>%
    summarize(count = n(), .groups = 'drop') %>%
    as.data.frame() %>% ungroup() %>%
    group_by(!!sym(cluster_res)) %>%
    mutate(perc =  round(count/sum(count)*100, 2)) #%>% ungroup()
head(cluster_summary)    
write.csv(cluster_summary, file.path(figs_dir_path_2, paste0("perc_samplesIncluster_summary_",cluster_res,".csv")))

p <- cluster_summary %>%
    ggplot(aes(x=.data[[cluster_res]], y=perc, fill=sample)) + 
    labs(x = "", y = "Cells (%)", fill="Sample") +
    geom_bar(stat="identity") +
    theme_classic() + 
    scale_fill_manual(values = sample_palette)
    # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave(file.path(figs_dir_path_2, paste0("barplot_clusterbysample_",cluster_res,".tiff")),
    plot = p, units="in", width=size*1.1, height=size*0.7, dpi=300, compression = 'lzw')  


#### =========================================== ####
#### UMAP of celltypes ####
#### =========================================== ####

makeUMAPPlot <- function(obj, ident, labels=FALSE) {
  DefaultAssay(obj) <- "integrated"
  Idents(obj) <- ident
  p <- DimPlot(obj, reduction = "umap", group.by = ident, label=labels, label.size = 3.5) + scale_colour_viridis_d()
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

p <- makeUMAPPlot(obj = seurat_obj_ge, ident = cluster_col)
ggsave(file.path(figs_dir_path, paste0("celltypes_DimPlot_nolabels",cluster_col,".tiff")),
  plot = p, units="in", width=size*1.6, height=size*1.6, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj = seurat_obj_ge, ident = cluster_col, labels = TRUE)
ggsave(file.path(figs_dir_path, paste0("celltypes_DimPlot",cluster_col,".tiff")),
  plot = p, units="in", width=size*1.6, height=size*1.6, dpi=300, compression = 'lzw')

p <- makeUMAPPlot(obj = seurat_obj_gte, ident = cluster_col)
ggsave(file.path(figs_dir_path_2, paste0("celltypes_DimPlot_nolabel",cluster_col,".tiff")),
  plot = p, units="in", width=size*1.6, height=size*1.6, dpi=300, compression = 'lzw')
p <- makeUMAPPlot(obj = seurat_obj_gte, ident = cluster_col, labels = TRUE)
ggsave(file.path(figs_dir_path_2, paste0("celltypes_DimPlot",cluster_col,".tiff")),
  plot = p, units="in", width=size*1.6, height=size*1.6, dpi=300, compression = 'lzw')

#### =========================================== ####
#### UMAP of known markers ####
#### =========================================== ####

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

# Create UMAPs per unique marker 
known_markers_unique <- unique(known_markers)
condition <- known_markers_unique %in% rownames(seurat_obj_ge)
for (gene in known_markers_unique[condition]) {
    print(gene)
    p <- makeUMAPPlot(seurat_obj_ge, gene)
    ggsave(file.path(figs_dir_path, paste0(gene, "_FeaturePlot.tiff")),
      plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
}

condition <- known_markers_unique %in% rownames(seurat_obj_gte)
for (gene in known_markers_unique[condition]) {
    print(gene)
    p <- makeUMAPPlot(seurat_obj_gte, gene)
    ggsave(file.path(figs_dir_path_2, paste0(gene, "_FeaturePlot.tiff")),
      plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
}

#### =========================================== ####
#### GSC & oRG Distribution if GBM dataset ####
#### =========================================== ####
if (grepl("gbm", subdir, fixed = TRUE)) {

  ## Export table counting cells per sample in each celltype + Bar plot with cell counts per cell type grouped by sample
  cluster_summary <- seurat_obj_ge@meta.data %>%
      group_by(sample, !!sym(cluster_col_gsc)) %>%
      summarize(count = n(), .groups = 'drop') %>%
      as.data.frame() %>% ungroup() %>%
      group_by(!!sym(cluster_col_gsc)) %>%
      mutate(perc =  round(count/sum(count)*100, 2)) #%>% ungroup()
  head(cluster_summary)    
  write.csv(cluster_summary, file.path(figs_dir_path, paste0("perc_samplesInCelltype_summary_",cluster_col_gsc,".csv")))

  p <- cluster_summary %>%
      ggplot(aes(x=!!sym(cluster_col_gsc), y=perc, fill=sample)) + 
      labs(x = "", y = "Cells (%)", fill="Sample") +
      geom_bar(stat="identity") +
      theme_classic() + 
      scale_fill_manual(values = sample_palette)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  ggsave(file.path(figs_dir_path, paste0("barplot_celltypeBysample",cluster_col_gsc,".tiff")),
      plot = p, units="in", width=size*1.3, height=size*0.9, dpi=300, compression = 'lzw')  


  cluster_summary <- seurat_obj_gte@meta.data %>%
      group_by(sample, !!sym(cluster_col_gsc)) %>%
      summarize(count = n(), .groups = 'drop') %>%
      as.data.frame() %>% ungroup() %>%
      group_by(!!sym(cluster_col_gsc)) %>%
      mutate(perc =  round(count/sum(count)*100, 2)) #%>% ungroup()
  head(cluster_summary)    
  write.csv(cluster_summary, file.path(figs_dir_path_2, paste0("perc_samplesInCelltype_summary_",cluster_col_gsc,".csv")))

  p <- cluster_summary %>%
      ggplot(aes(x=.data[[cluster_col_gsc]], y=perc, fill=sample)) + 
      labs(x = "", y = "Cells (%)", fill="Sample") +
      geom_bar(stat="identity") +
      theme_classic() + 
      scale_fill_manual(values = sample_palette)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
  ggsave(file.path(figs_dir_path_2, paste0("barplot_celltypeBysample_",cluster_col_gsc,".tiff")),
      plot = p, units="in", width=size*1.3, height=size*0.8, dpi=300, compression = 'lzw')  

  # Dotplots (expression)
  DefaultAssay(seurat_obj_ge) <- "RNA"
  Idents(seurat_obj_ge) <- cluster_col_gsc

  condition <- known_markers %in% rownames(seurat_obj_ge)
  print(paste("Length of unique features to create DotPlot", length(unique(known_markers))))

  p <- DotPlot(seurat_obj_ge, features = rev(unique(known_markers[condition]))) + scale_colour_gradient2() + 
        xlab("") + ylab("Cell Types") + coord_flip() + 
        theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
  ggsave(file.path(figs_dir_path, paste0("celltypeXuniqueMarkers_DotPlot_",cluster_col_gsc,".tiff")),
      plot = p, units="in", width=size*1.8, height=size*1.8, dpi=300, compression = 'lzw')

  DefaultAssay(seurat_obj_gte) <- "RNA"
  Idents(seurat_obj_gte) <- cluster_col_gsc

  condition <- known_markers %in% rownames(seurat_obj_gte)
  print(paste("Length of unique features to create DotPlot", length(unique(known_markers))))

  p <- DotPlot(seurat_obj_gte, features = rev(unique(known_markers[condition]))) + scale_colour_gradient2() + 
        xlab("") + ylab("Cell Types") + coord_flip() + 
        theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
  ggsave(file.path(figs_dir_path_2, paste0("celltypeXuniqueMarkers_DotPlot_",cluster_col_gsc,".tiff")),
      plot = p, units="in", width=size*1.8, height=size*1.8, dpi=300, compression = 'lzw')

  # UMAPs of Gene Markers
  condition <- s.genes %in% rownames(seurat_obj_ge)
  for (gene in s.genes[condition]) {
      print(gene)
      p <- makeUMAPPlot(seurat_obj_ge, gene)
      ggsave(file.path(figs_dir_path, paste0(gene, "_CC.S_FeaturePlot.tiff")),
        plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
  }

  condition <- g2m.genes %in% rownames(seurat_obj_ge)
  for (gene in g2m.genes[condition]) {
      print(gene)
      p <- makeUMAPPlot(seurat_obj_ge, gene)
      ggsave(file.path(figs_dir_path, paste0(gene, "_CC.G2M_FeaturePlot.tiff")),
        plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
  }
  gsc_genes <- c('PROM1','FUT4','L1CAM')
  condition <- gsc_genes %in% rownames(seurat_obj_ge)
  for (gene in gsc_genes[condition]) {
      print(gene)
      p <- makeUMAPPlot(seurat_obj_ge, gene)
      ggsave(file.path(figs_dir_path, paste0(gene, "_GSC_FeaturePlot.tiff")),
        plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
  }

  condition <- s.genes %in% rownames(seurat_obj_gte)
  for (gene in s.genes[condition]) {
      print(gene)
      p <- makeUMAPPlot(seurat_obj_gte, gene)
      ggsave(file.path(figs_dir_path_2, paste0(gene, "_CC.S_FeaturePlot.tiff")),
        plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
  }
  condition <- g2m.genes %in% rownames(seurat_obj_gte)
  for (gene in g2m.genes[condition]) {
      print(gene)
      p <- makeUMAPPlot(seurat_obj_gte, gene)
      ggsave(file.path(figs_dir_path_2, paste0(gene, "_CC.G2M_FeaturePlot.tiff")),
        plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
  }
  condition <- gsc_genes %in% rownames(seurat_obj_gte)
  for (gene in gsc_genes[condition]) {
      print(gene)
      p <- makeUMAPPlot(seurat_obj_gte, gsc_genes)
      ggsave(file.path(figs_dir_path_2, paste0(gene, "_GSC_FeaturePlot.tiff")),
        plot = p, units="in", width=size*0.8, height=size*0.8, dpi=300, compression = 'lzw')
  }

  # UMAP of GSC and oRGs

  p <- DimPlot(seurat_obj_ge, reduction = "umap", 
              cells.highlight = WhichCells(object = seurat_obj_ge, cells = barcode_org) )
  ggsave(file.path(figs_dir_path, paste0("org_DimPlot.tiff")),
    plot = p, units="in", width=size*2.0, height=size*1.6, dpi=300, compression = 'lzw')
  p <- DimPlot(seurat_obj_gte, reduction = "umap", 
              cells.highlight = WhichCells(object = seurat_obj_gte, cells = barcode_org) )
  ggsave(file.path(figs_dir_path_2, paste0("org_DimPlot.tiff")),
    plot = p, units="in", width=size*2.0, height=size*1.6, dpi=300, compression = 'lzw')

  p <- DimPlot(seurat_obj_ge, reduction = "umap", 
              cells.highlight = WhichCells(object = seurat_obj_ge, cells = barcode_gsc) )
  ggsave(file.path(figs_dir_path, paste0("gsc_DimPlot.tiff")),
    plot = p, units="in", width=size*2.0, height=size*1.6, dpi=300, compression = 'lzw')
  p <- DimPlot(seurat_obj_gte, reduction = "umap", 
              cells.highlight = WhichCells(object = seurat_obj_gte, cells = barcode_gsc) )
  ggsave(file.path(figs_dir_path_2, paste0("gsc_DimPlot.tiff")),
    plot = p, units="in", width=size*2.0, height=size*1.6, dpi=300, compression = 'lzw')

  #### =========================================== ####
  #### Stem Cell Scores ####
  #### =========================================== ####

  table(seurat_obj_ge@meta.data[c('Phase')])
  #    G1   G2M     S
  # 21983  4010  6090

  table(seurat_obj_ge@meta.data[c('Phase','gsc')])
  #      gsc
  # Phase     0     1
  #   G1  20239  1744
  #   G2M  3755   255
  #   S    5469   621

  table(seurat_obj_ge@meta.data[c('Phase','oRG')])
  #      oRG
  # Phase     0     1
  #   G1  21685   298
  #   G2M  3990    20
  #   S    5995    95

  table(seurat_obj_ge@meta.data[c(cluster_res,'oRG')])
  #                       oRG
  # integrated_snn_res.0.6    0    1
  #                     0  4344   57
  #                     1  3945   84
  #                     2  3552   98
  #                     3  3035   98
  #                     4  2293    1
  #                     5  2213    0
  #                     6  1806    1
  #                     7  1569   27
  #                     7  1569   27
  #                     8  1424    0
  #                     9  1331    2
  #                     10 1246   32
  #                     11  950    6
  #                     12  934    0
  #                     13  896    0
  #                     14  842    2
  #                     15  639    0
  #                     16  397    5
  #                     17  184    0
  #                     18   70    0

  calcProportions <- function (meta, col1, col2) {
    grouped_counts_table <- table(meta[c(col1, col2)])
    grouped_counts_table <- cbind(grouped_counts_table, prop.table(grouped_counts_table, margin=1)*100)
    return(grouped_counts_table)
  }

  stats_df <- seurat_obj_ge@meta.data %>% 
    group_by(!!sym(cluster_res)) %>% 
    summarise(mean.s = mean(S.Score), sd.s = sd(S.Score, na.rm = TRUE), 
              mean.g2m = mean(G2M.Score), sd.g2m = sd(G2M.Score, na.rm = TRUE),
              t_pval = t.test(S.Score,G2M.Score, var.equal=TRUE)$p.value, #two-sided default
              t_statistic = t.test(S.Score,G2M.Score, var.equal=TRUE)$statistic,
              t_conf.int_min = t.test(S.Score,G2M.Score, var.equal=TRUE)$conf.int[1], # conf.level 0.95
              t_conf.int_max = t.test(S.Score,G2M.Score, var.equal=TRUE)$conf.int[2], # conf.level 0.95
              welch_pval = t.test(S.Score,G2M.Score, var.equal=FALSE)$p.value, #two-sided default
              welch_statistic = t.test(S.Score,G2M.Score, var.equal=FALSE)$statistic,
              welch_conf.int_min = t.test(S.Score,G2M.Score, var.equal=FALSE)$conf.int[1], # conf.level 0.95
              welch_conf.int_max = t.test(S.Score,G2M.Score, var.equal=FALSE)$conf.int[2], # conf.level 0.95
              .groups = "drop") #two-sided default

  var_tests <- seurat_obj_ge@meta.data %>% select(!!sym(cluster_res), S.Score, G2M.Score) %>%
    melt(value.name = "scores", id = cluster_res) %>%
    group_by(!!sym(cluster_res)) %>% 
    summarise( levene_pval = leveneTest(scores~variable, center=mean)$`Pr(>F)`[1], .groups = "drop" )

  grouped_counts_table <- calcProportions(seurat_obj_ge@meta.data, cluster_res,'Phase')
  summary_df <- as.data.frame.matrix(grouped_counts_table)
  colnames(summary_df) <- c("ncells.G1", "ncells.G2M", "ncells.S", "perc.G1", "perc.G2M", "perc.S")
  write.csv(cbind(stats_df, var_tests[,2], summary_df),file.path(figs_dir_path,paste0("cellcyclescore_summarystats",cluster_res,".csv")), row.names = FALSE)
  write.csv(cbind(stats_df, var_tests[,2], summary_df),file.path(figs_dir_path_2,paste0("cellcyclescore_summarystats",cluster_res,".csv")), row.names = FALSE)

  print("Exporting summary csvs for gsc and org labeled cells GE")
  write.csv(calcProportions(seurat_obj_ge@meta.data, cluster_res, 'gsc'), file.path(figs_dir_path, paste0("summary_gsc_res0.6VSgsc.csv")))
  write.csv(calcProportions(seurat_obj_ge@meta.data, cluster_res,'oRG'), file.path(figs_dir_path, paste0("summary_gsc_res0.6VSoRG.csv")))
  write.csv(calcProportions(seurat_obj_ge@meta.data, cluster_res,'Phase'), file.path(figs_dir_path, paste0("summary_gsc_res0.6VSccphase.csv")))
  write.csv(calcProportions(seurat_obj_ge@meta.data, cluster_res,cluster_col_gsc), file.path(figs_dir_path, paste0("summary_gsc_res0.6VS",cluster_col_gsc,".csv")))
  write.csv(calcProportions(seurat_obj_ge@meta.data, 'Phase','gsc'), file.path(figs_dir_path, paste0("summary_gsc_ccphaseVSgsc.csv")))
  write.csv(calcProportions(seurat_obj_ge@meta.data, 'Phase','oRG'), file.path(figs_dir_path, paste0("summary_gsc_ccphaseVSoRG.csv")))
  write.csv(calcProportions(seurat_obj_ge@meta.data, 'sample','oRG'), file.path(figs_dir_path, paste0("summary_gsc_ccphaseVSoRG.csv")))
  write.csv(calcProportions(seurat_obj_ge@meta.data, 'sample',cluster_col_gsc), file.path(figs_dir_path, paste0("summary_gsc_sampleVS",cluster_col_gsc,".csv")))

  print("Exporting summary csvs for gsc and org labeled cells GTE")
  write.csv(calcProportions(seurat_obj_gte@meta.data, cluster_res, 'gsc'), file.path(figs_dir_path_2, paste0("summary_gsc_res0.6VSgsc.csv")))
  write.csv(calcProportions(seurat_obj_gte@meta.data, cluster_res,'oRG'), file.path(figs_dir_path_2, paste0("summary_gsc_res0.6VSoRG.csv")))
  write.csv(calcProportions(seurat_obj_gte@meta.data, cluster_res,'Phase'), file.path(figs_dir_path_2, paste0("summary_gsc_res0.6VSccphase.csv")))
  write.csv(calcProportions(seurat_obj_gte@meta.data, cluster_res,cluster_col_gsc), file.path(figs_dir_path_2, paste0("summary_gsc_res0.6VS",cluster_col_gsc,".csv")))
  write.csv(calcProportions(seurat_obj_gte@meta.data, 'Phase','gsc'), file.path(figs_dir_path_2, paste0("summary_gsc_ccphaseVSgsc.csv")))
  write.csv(calcProportions(seurat_obj_gte@meta.data, 'Phase','oRG'), file.path(figs_dir_path_2, paste0("summary_gsc_ccphaseVSoRG.csv")))
  write.csv(calcProportions(seurat_obj_gte@meta.data, 'sample','oRG'), file.path(figs_dir_path_2, paste0("summary_gsc_ccphaseVSoRG.csv")))
  write.csv(calcProportions(seurat_obj_gte@meta.data, 'sample',cluster_col_gsc), file.path(figs_dir_path_2, paste0("summary_gsc_sampleVS",cluster_col_gsc,".csv")))

  # Create Barplot for CellCycleScore (grouped by cluster)
  library(reshape2)
  print(stats_df)
  if (colnames(stats_df)[1] != cluster_res) { # convert row names to columns if column name doesn't exist
    stats_df <- rownames_to_column(stats_df, var = cluster_res) 
  }
  summary_means <- stats_df %>% select(!!sym(cluster_res), mean.s, mean.g2m) %>%
    melt(value.name = "means", id = cluster_res)
  summary_sd <- stats_df %>% select(!!sym(cluster_res), sd.s, sd.g2m) %>%
    melt(value.name = "sds", id = cluster_res)
  summary_stats <- cbind(summary_means, summary_sd[,3])
  colnames(summary_stats)[2] <- "phase"
  colnames(summary_stats)[4] <- "sds"

  p <- summary_stats %>% ggplot(aes(x=.data[[cluster_res]], y=means, fill=factor(phase)) ) +  
    geom_bar(position=position_dodge(), stat="identity", colour='black') +
    scale_fill_viridis_d(labels=c('Mean(S Phase Score)', 'Mean(G2/M Phase Score)')) + theme_classic()+
    labs(x="", y = "Mean(Cell Cycle Phase Score)", fill="Phase") + #NoLegend() + 
    guides(fill = guide_legend(position = "inside")) +
    theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
          axis.title.y = element_text(size = 14),
          legend.position.inside = c(0.3, 0.9)) +
    geom_errorbar(aes(ymin=means-sds, ymax=means+sds), width=.2, position = position_dodge(0.9)) 
  ggsave(file.path(figs_dir_path, paste0("barplot_ccphasescore_",cluster_res,".tiff")),
      plot = p, units="in", width=size*1, height=size*0.7, dpi=300, compression = 'lzw')  

  # Create Barplot for CellCycleScore (grouped by gsctypes)
  stats_df <- seurat_obj_ge@meta.data %>% 
    group_by(!!sym(cluster_col_gsc)) %>% 
    summarise(mean.s = mean(S.Score), sd.s = sd(S.Score, na.rm = TRUE), 
              mean.g2m = mean(G2M.Score), sd.g2m = sd(G2M.Score, na.rm = TRUE),
              t_pval = t.test(S.Score,G2M.Score, var.equal=TRUE)$p.value, #two-sided default
              t_statistic = t.test(S.Score,G2M.Score, var.equal=TRUE)$statistic,
              t_conf.int_min = t.test(S.Score,G2M.Score, var.equal=TRUE)$conf.int[1], # conf.level 0.95
              t_conf.int_max = t.test(S.Score,G2M.Score, var.equal=TRUE)$conf.int[2], # conf.level 0.95
              welch_pval = t.test(S.Score,G2M.Score, var.equal=FALSE)$p.value, #two-sided default
              welch_statistic = t.test(S.Score,G2M.Score, var.equal=FALSE)$statistic,
              welch_conf.int_min = t.test(S.Score,G2M.Score, var.equal=FALSE)$conf.int[1], # conf.level 0.95
              welch_conf.int_max = t.test(S.Score,G2M.Score, var.equal=FALSE)$conf.int[2], # conf.level 0.95
              .groups = "drop") #two-sided default

  var_tests <- seurat_obj_ge@meta.data %>% select(!!sym(cluster_col_gsc), S.Score, G2M.Score) %>%
    melt(value.name = "scores", id = cluster_col_gsc) %>%
    group_by(!!sym(cluster_col_gsc)) %>% 
    summarise( levene_pval = leveneTest(scores~variable, center=mean)$`Pr(>F)`[1], .groups = "drop" )

  grouped_counts_table <- calcProportions(seurat_obj_ge@meta.data, cluster_col_gsc,'Phase')
  summary_df <- as.data.frame.matrix(grouped_counts_table)
  colnames(summary_df) <- c("ncells.G1", "ncells.G2M", "ncells.S", "perc.G1", "perc.G2M", "perc.S")
  write.csv(cbind(stats_df, var_tests[,2], summary_df),file.path(figs_dir_path,paste0("cellcyclescore_summarystats_",cluster_col_gsc,".csv")), row.names = FALSE)
  write.csv(cbind(stats_df, var_tests[,2], summary_df),file.path(figs_dir_path_2,paste0("cellcyclescore_summarystats_",cluster_col_gsc,".csv")), row.names = FALSE)

  library(reshape2)
  print(stats_df)
  if (colnames(stats_df)[1] != cluster_col_gsc) { # convert row names to columns if column name doesn't exist
    stats_df <- rownames_to_column(stats_df, var = cluster_col_gsc) 
  }
  summary_means <- stats_df %>% select(!!sym(cluster_col_gsc), mean.s, mean.g2m) %>%
    melt(value.name = "means", id = cluster_col_gsc)
  summary_sd <- stats_df %>% select(!!sym(cluster_col_gsc), sd.s, sd.g2m) %>%
    melt(value.name = "sds", id = cluster_col_gsc)
  summary_stats <- cbind(summary_means, summary_sd[,3])
  colnames(summary_stats)[2] <- "phase"
  colnames(summary_stats)[4] <- "sds"

  p <- summary_stats %>% ggplot(aes(x=.data[[cluster_col_gsc]], y=means, fill=factor(phase)) ) +  
    geom_bar(position=position_dodge(0.9), stat="identity", colour='black') +
    scale_fill_viridis_d(labels=c('Mean(S Phase Score)', 'Mean(G2/M Phase Score)')) + theme_classic()+
    labs(x="", y = "Mean(Cell Cycle Phase Score)", fill="Phase") + #NoLegend() + 
    guides(fill = guide_legend(position = "inside")) +
    theme(axis.text.x = element_text(size=12, angle = 40, vjust = 1, hjust=1), 
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size = 14),
          legend.position.inside = c(0.5, 0.8)) +
    geom_errorbar(aes(ymin=means-sds, ymax=means+sds), width=.2, position = position_dodge(0.9)) 
  ggsave(file.path(figs_dir_path, paste0("barplot_ccphasescore_",cluster_col_gsc,".tiff")),
      plot = p, units="in", width=size*1.3, height=size*0.9, dpi=300, compression = 'lzw')  
  ggsave(file.path(figs_dir_path_2, paste0("barplot_ccphasescore_",cluster_col_gsc,".tiff")),
      plot = p, units="in", width=size*1.3, height=size*0.9, dpi=300, compression = 'lzw')  

  #### =========================================== ####
  #### Stem Cell Scores (cc celltypes) ####
  #### =========================================== ####

  # Create Barplot for CellCycleScore (grouped by gsctypes)
  stats_df <- seurat_obj_ge@meta.data %>% 
    group_by(!!sym(cluster_col_gsc_cc)) %>% 
    summarise(mean.s = mean(S.Score), sd.s = sd(S.Score, na.rm = TRUE), 
              mean.g2m = mean(G2M.Score), sd.g2m = sd(G2M.Score, na.rm = TRUE),
              t_pval = t.test(S.Score,G2M.Score, var.equal=TRUE)$p.value, #two-sided default
              t_statistic = t.test(S.Score,G2M.Score, var.equal=TRUE)$statistic,
              t_conf.int_min = t.test(S.Score,G2M.Score, var.equal=TRUE)$conf.int[1], # conf.level 0.95
              t_conf.int_max = t.test(S.Score,G2M.Score, var.equal=TRUE)$conf.int[2], # conf.level 0.95
              welch_pval = t.test(S.Score,G2M.Score, var.equal=FALSE)$p.value, #two-sided default
              welch_statistic = t.test(S.Score,G2M.Score, var.equal=FALSE)$statistic,
              welch_conf.int_min = t.test(S.Score,G2M.Score, var.equal=FALSE)$conf.int[1], # conf.level 0.95
              welch_conf.int_max = t.test(S.Score,G2M.Score, var.equal=FALSE)$conf.int[2], # conf.level 0.95
              .groups = "drop") #two-sided default

  var_tests <- seurat_obj_ge@meta.data %>% select(!!sym(cluster_col_gsc_cc), S.Score, G2M.Score) %>%
    melt(value.name = "scores", id = cluster_col_gsc_cc) %>%
    group_by(!!sym(cluster_col_gsc_cc)) %>% 
    summarise( levene_pval = leveneTest(scores~variable, center=mean)$`Pr(>F)`[1], .groups = "drop" )

  grouped_counts_table <- calcProportions(seurat_obj_ge@meta.data, cluster_col_gsc_cc,'Phase')
  summary_df <- as.data.frame.matrix(grouped_counts_table)
  colnames(summary_df) <- c("ncells.G1", "ncells.G2M", "ncells.S", "perc.G1", "perc.G2M", "perc.S")
  write.csv(cbind(stats_df, var_tests[,2], summary_df),file.path(figs_dir_path,paste0("cellcyclescore_summarystats_",cluster_col_gsc_cc,".csv")), row.names = FALSE)
  write.csv(cbind(stats_df, var_tests[,2], summary_df),file.path(figs_dir_path_2,paste0("cellcyclescore_summarystats_",cluster_col_gsc_cc,".csv")), row.names = FALSE)

  library(reshape2)
  print(stats_df)
  if (colnames(stats_df)[1] != cluster_col_gsc_cc) { # convert row names to columns if column name doesn't exist
    stats_df <- rownames_to_column(stats_df, var = cluster_col_gsc_cc) 
  }
  summary_means <- stats_df %>% select(!!sym(cluster_col_gsc_cc), mean.s, mean.g2m) %>%
    melt(value.name = "means", id = cluster_col_gsc_cc)
  summary_sd <- stats_df %>% select(!!sym(cluster_col_gsc_cc), sd.s, sd.g2m) %>%
    melt(value.name = "sds", id = cluster_col_gsc_cc)
  summary_stats <- cbind(summary_means, summary_sd[,3])
  colnames(summary_stats)[2] <- "phase"
  colnames(summary_stats)[4] <- "sds"

  p <- summary_stats %>% ggplot(aes(x=.data[[cluster_col_gsc_cc]], y=means, fill=factor(phase)) ) +  
    geom_bar(position=position_dodge(0.9), stat="identity", colour='black') +
    scale_fill_viridis_d(labels=c('Mean(S Phase Score)', 'Mean(G2/M Phase Score)')) + theme_classic()+
    labs(x="", y = "Mean(Cell Cycle Phase Score)", fill="Phase") + #NoLegend() + 
    guides(fill = guide_legend(position = "inside")) +
    theme(axis.text.x = element_text(size=12, angle = 40, vjust = 1, hjust=1), 
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size = 14),
          legend.position.inside = c(0.5, 0.8)) +
    geom_errorbar(aes(ymin=means-sds, ymax=means+sds), width=.2, position = position_dodge(0.9)) 
  ggsave(file.path(figs_dir_path, paste0("barplot_ccphasescore_",cluster_col_gsc_cc,".tiff")),
      plot = p, units="in", width=size*1.3, height=size*0.9, dpi=300, compression = 'lzw')  
  ggsave(file.path(figs_dir_path_2, paste0("barplot_ccphasescore_",cluster_col_gsc_cc,".tiff")),
      plot = p, units="in", width=size*1.3, height=size*0.9, dpi=300, compression = 'lzw')  
}

### End of Script
sessionInfo()