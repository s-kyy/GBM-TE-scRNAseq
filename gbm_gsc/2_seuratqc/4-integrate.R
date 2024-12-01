#!usr/bin/env Rscript

#### =========================================== ####
#### Verify Args ####
#### =========================================== ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("At least 2 filepaths names must be supplied: [dataset/ge.rds] [cellcycle_genes.csv]. Optionally include path to an output folder [output_path]", call.=FALSE)
} else if (length(args)>=2) {
  # verify filepaths
  if (file.exists(args[1]) & file.exists(args[2])){ 
    obj_path <- args[1] 
    ccgenes_path <- args[2]
    filename <- basename(path_ext_remove(obj_path))
    parent_dir_path <- dirname(obj_path)
    parent_dir_name <- basename(parent_dir_path)
  } else {
    stop("one or more filepaths do not exist. Closing script...", call=FALSE)
  }
  # Optional arguements
  if (length(args) == 3 & dir.exists(args[3])) {
    parent_dir_path <- args[3]
    parent_dir_name <- basename(parent_dir_path)
  } else if (length(args) == 3 & !dir.exists(args[3])) {
    print("Output directory does not exist, creating new output directory...")
    dir.create(parent_dir_path, recursive=T)
  }
} 

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/tcga-gbm-R4/renv/library/R-4.0/x86_64-pc-linux-gnu", .libPaths()))
.libPaths()

set.seed(108)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 

set.seed(108)

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
obj <- readRDS(obj_path) 
cc.genes <- read.csv(ccgenes_path)
print(paste("First ten genes out of", length(cc.genes$s.genes),
            "associated to S phase in mitosis", 
            paste(head(cc.genes$s.genes,10))))
print(paste("First ten genes out of", length(cc.genes$g2m.genes),
            "associated to G2/M phase in mitosis", 
            paste(head(cc.genes$g2m.genes,10))))
print(paste(length(intersect(cc.genes$s.genes, cc.genes$g2m.genes)),
            "genes associated to both S phase and G2/M phase in mitosis",
            paste(intersect(cc.genes$s.genes, cc.genes$g2m.genes))))
ndims <- 25
size <- 5

figs_dir_path <- file.path(parent_dir_path, "figs")

ifelse(!dir.exists(figs_dir_path),
        dir.create(figs_dir_path,recursive=T),
        "Directory Exists")

#### =========================================== ####
#### Normalize & Find Variable Features each dataset individually ####
#### =========================================== ####

# Correct for inter-sample variability 
# (i.e. mainly from different tumour donors) 
print(paste("unique samples to split by:", unique(obj$sample)))
obj_list <- SplitObject(object=obj, split.by="sample")
print(paste("Length of list", length(obj_list)))
rm("obj")
gc()

# Normalize the RNA reads per tumor sample via Log Normalization method
p_list <- vector(mode="list",length(obj_list))
names(p_list) <- sort(names(obj_list))

for(i in 1:length(obj_list)){
    print(paste0("Processing Sample: ", obj_list[[i]]$sample[1]))
    DefaultAssay(obj_list[[i]]) <- "RNA"
    # default: LogNormalize, save to : $RNA@data
    obj_list[[i]] <- NormalizeData(obj_list[[i]], verbose = TRUE) 
    obj_list[[i]] <- FindVariableFeatures(obj_list[[i]], 
                              selection.method="vst",
                              nfeatures=2500, #default=2000
                              verbose = TRUE)
    #Create Variable Feature Plots
    top10 <- head(VariableFeatures(obj_list[[i]]), 10)
    p <- VariableFeaturePlot(obj_list[[i]])
    p <- LabelPoints(plot = p, points = top10, repel = TRUE)
    ggsave(file.path(figs_dir_path, paste0(filename,obj_list[[i]]$sample[1],"_varfeatgenes.tiff")), 
          plot = p, units="in", width=size*1.1, height=size*0.8, dpi=300, compression = 'lzw')
    p_list[[obj_list[[i]]$sample[1]]] <- p
}

saveRDS(p_list, file=file.path(figs_dir_path,paste0(filename,"_allvarfeatplots.rds")))

print("Variable Features computed")

#### =========================================== ####
#### Integrate Seurat Objects together: 30min-1h30min ####
#### =========================================== ####
## Anchors across all samples will be used to integrated each sample
## This step helps reduce variation introduced across donor samples 
## (i.e. variation in tumour extraction, experimental processing of samples, 
## gender-based differences and so on)
## Stuck with CCA method to find anchors, because the datasets are a manageable size (<50K cells)
obj_anchors <- FindIntegrationAnchors(
                  object.list = obj_list, 
                  dims=1:ndims,
                  normalization.method = "LogNormalize", 
                  anchor.features = 2500, #match number in FindVariableFeatures
                  )
obj_integrated <- IntegrateData(anchorset = obj_anchors, dims=1:ndims)
print("Integration Complete")

saveRDS(obj_integrated, file = file.path(parent_dir_path, paste0(filename, "_integrated.rds")) )
print(paste("Saved integrated seurat objects to",parent_dir_path))

rm("obj_list")
rm("obj_anchors")
gc()

#### =========================================== ####
#### Scale, PCA, UMAP ~5min, <200Mb ####
#### =========================================== ####
## save to $RNA@scale.data
DefaultAssay(obj_integrated) <- "RNA"

## ========================================= ##
##  Regress out cell cycle effects completely ##
## ========================================= ##
# Following Satija Lab Vignette:  https://satijalab.org/seurat/archive/v4.3/cell_cycle_vignette
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
# other markers: https://www.nature.com/articles/s41698-022-00302-7
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

obj_integrated <- CellCycleScoring(obj_integrated, 
                                    s.features = s.genes, 
                                    g2m.features = g2m.genes, 
                                    set.ident = FALSE)
# Visualize the distribution of cell cycle markers across
p <- RidgePlot(obj_integrated, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
ggsave(file.path(figs_dir_path, paste0(filename,"_cellcycleRidgePlot.tiff")), 
      plot = p, units="in", width=size*1.5, height=size*1, dpi=300, compression = 'lzw')

# Regress out cell cycle score
DefaultAssay(obj_integrated) <- "integrated"
obj_integrated_cc <- ScaleData(obj_integrated, 
                            vars.to.regress = c("S.Score", "G2M.Score"), 
                            features = rownames(obj_integrated))

# Run PCA and UMAP ~20min
obj_integrated_cc <- RunPCA(obj_integrated_cc, 
                          features = VariableFeatures(obj_integrated_cc),
                          npcs = ndims, 
                          verbose = FALSE)
obj_integrated_cc <- RunUMAP(obj_integrated_cc,
                          reduction = "pca", 
                          dims = 1:ndims, 
                          umap.method = "uwot", 
                          metric = "cosine")

saveRDS(obj_integrated_cc, file = file.path(parent_dir_path, paste0(filename, "_integrated_regressCC.rds")) )
print(paste("ScaleData, RunPCA, and RunUMAP complete.\n Saved Seurat objects to",parent_dir_path))
rm("obj_integrated_cc")

## ========================================= ##
## Regress out differences in cell cycle effects ##
## ========================================= ##
# Regress out differences within cycling cells and retain features 
# differentiating proliferating cells from non-proliferating cells
DefaultAssay(obj_integrated) <- "integrated"
obj_integrated$CC.Difference <- obj_integrated$S.Score - obj_integrated$G2M.Score
obj_integrated_cc <- ScaleData(obj_integrated, 
                            vars.to.regress = "CC.Difference", 
                            features = rownames(obj_integrated), 
                            verbose = FALSE)

obj_integrated_cc@meta.data$sample <- factor(obj_integrated_cc@meta.data$sample, 
                                          levels = sort(unique(obj_integrated_cc@meta.data$sample)))
# Run PCA and UMAP ~20min
obj_integrated_cc <- RunPCA(obj_integrated_cc, 
                          features = VariableFeatures(obj_integrated_cc),
                          npcs = ndims, 
                          verbose = FALSE)

obj_integrated_cc <- RunUMAP(obj_integrated_cc, 
                          reduction = "pca", 
                          dims = 1:ndims, 
                          umap.method = "uwot", 
                          metric = "cosine")

saveRDS(obj_integrated_cc, file = file.path(parent_dir_path, paste0(filename, "_integrated_regressCCdiff.rds")) )
print(paste("ScaleData, RunPCA, and RunUMAP complete.\n Saved Seurat objects to",parent_dir_path))

## ========================================= ##
## Scale and Reduce Dimensionality ##
## ========================================= ##
# obj_integrated <- ScaleData(object = obj_integrated, verbose = FALSE) 

# obj_integrated@meta.data$sample <- factor(obj_integrated@meta.data$sample, 
#                                           levels = sort(unique(obj_integrated@meta.data$sample)))
# # Run PCA and UMAP ~20min
# obj_integrated <- RunPCA(obj_integrated, 
#                           npcs = ndims, 
#                           verbose = FALSE)

# obj_integrated <- RunUMAP(obj_integrated, 
#                           reduction = "pca", 
#                           dims = 1:ndims, 
#                           umap.method = "uwot", 
#                           metric = "cosine")
# saveRDS(obj_integrated, file = file.path(parent_dir_path, paste0(filename, "_integrated_regressCCdiff.rds")) )
# print(paste("ScaleData, RunPCA, and RunUMAP complete.\n Saved Seurat objects to",parent_dir_path))

#### End of Script ####
sessionInfo()
