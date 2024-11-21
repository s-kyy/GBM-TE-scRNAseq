#!usr/bin/env Rscript

#### =========================================== ####
#### Verify Args ####
#### =========================================== ####
library(fs) #path manipulation

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<=1) {
  stop("At least 2 filepaths names must be supplied: [dataset1/ge.rds] [dataset2/gte.rds]. Optionally include path to an output folder [output_path]", call.=FALSE)
} else if (length(args)>=2) {
  # verify filepaths
  if (file.exists(args[1]) & file.exists(args[2])){ 
    path_ge <- args[1]  
    path_gte <- args[2] 
    ge_filename <- basename(path_ext_remove(path_ge))
    gte_filename <- basename(path_ext_remove(path_gte))
    path_to_object <- dirname(path_ge)
    parent_dir_name <- basename(path_to_object)
  } else {
    stop("filepaths for dataset1 and/or dataset2 does not exist. Closing script...", call=FALSE)
  }
  # Optional arguements
  if (length(args) == 3 & dir.exists(args[3])) {
    path_to_object <- args[3]
    parent_dir_name <- basename(path_to_object)
  } else if (length(args) == 3 & !dir.exists(args[3])) {
    cat("Output directory does not exist, creating new output directory...")
    dir.create(path_to_object, recursive=T)
  }
} 

#### =========================================== ####
#### Import Packages ####
#### =========================================== ####
set.seed(34)

library(Seurat)
library(Matrix)
library(ggplot2)
library(tidyverse) 
library(scales) #plot axis manipulation

set.seed(34)

#### =========================================== ####
#### Load Datasets ####
#### =========================================== ####
ge <- readRDS(path_ge) 
gte <- readRDS(path_gte) 

#### =========================================== ####
#### Filter Low Quality Cells ####
#### =========================================== ####
filt.ge <- subset(x = ge, subset= (nUMI >= 1000) & 
                            (nGene >= 300) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

barcode <- colnames(filt.ge)
filt.gte <- subset(gte, cells=barcode)

cat(paste0("\n ge contains", 
           dim(filt.ge), " genes x cells\n"))
cat(paste0("Filtered out ", 
           dim(ge)[2] - dim(filt.ge)[2], " cells in ge\n" ))

cat(paste0("\n gte contains", 
           dim(filt.gte), " genes x cells\n"))
cat(paste0("Filtered out ", 
           dim(gte)[2] - dim(filt.gte)[2], " cells in ge\n" ))

saveRDS(filt.gte, file = file.path(path_to_object, paste0(gte_filename, "_qc.rds")) )
saveRDS(filt.ge,  file = file.path(path_to_object, paste0(ge_filename, "_qc.rds")) )
cat("Filtered & saved seurat objects to",path_to_object,"\n")

rm("ge")
rm("gte")
gc()

#### =========================================== ####
#### Normalize each dataset individually ####
#### =========================================== ####
gte.temp <- CreateAssayObject(counts = filt.gte@assays[["RNA"]]@counts)

gte.temp <- NormalizeData(gte.temp, 
                          verbose = T, 
                          normalization.method = "LogNormalize") %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData()

filt.gte@assays[["RNA"]]@counts <- gte.temp@counts
filt.gte@assays[["RNA"]]@scale.data <- gte.temp@scale.data

ge.temp  <- CreateAssayObject(counts = filt.ge@assays[["RNA"]]@counts)
ge.temp <- NormalizeData(ge.temp, 
                          verbose = T, 
                          normalization.method = "LogNormalize") %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData()

filt.ge@assays[["RNA"]]@counts <- ge.temp@counts
filt.ge@assays[["RNA"]]@scale.data <- ge.temp@scale.data

rm("ge.temp")
rm("gte.temp")
gc()

gte.list <- SplitObject(object=filt.gte, split.by="orig.ident")
ge.list <- SplitObject(object=filt.ge, split.by="orig.ident")

rm("filt.gte")
rm("filt.ge")
gc()


# Change active.ident as orig.ident for GBM samples
for (i in 1:length(gte.list)) {
    gte.list[[i]]@active.ident <- factor(gte.list[[i]]@meta.data$orig.ident)
    ge.list[[i]]@active.ident <- factor(ge.list[[i]]@meta.data$orig.ident)
}

# Normalize the RNA reads per tumor sample via Log Normalization method

for(i in 1:length(ge.list)){
    
    print(paste0("Processing Sample Number: ",i))
    
    DefaultAssay(gte.list[[i]]) <- "RNA"
    DefaultAssay(ge.list[[i]]) <- "RNA"
    
    gte.list[[i]] <- NormalizeData(gte.list[[i]],verbose=TRUE)
    ge.list[[i]]  <- NormalizeData(ge.list[[i]],verbose=TRUE)
    
    gte.list[[i]] <- FindVariableFeatures(gte.list[[i]],
                                          selection.method="vst",
                                          nfeatures=2500,
                                          verbose=TRUE)
    ge.list[[i]]  <- FindVariableFeatures(ge.list[[i]],
                                          selection.method="vst",
                                          nfeatures=2500,
                                          verbose=TRUE)
}
cat("Most Variable Features computed\n")

#### =========================================== ####
#### Integrate Seurat Objects together (memory: 60Gb) ~1h30min ####
#### =========================================== ####
gte.anchors <- FindIntegrationAnchors(gte.list,dims=1:20)
gte.integrated <- IntegrateData(anchorset = gte.anchor,dims=1:20)

ge.anchors <- FindIntegrationAnchors(ge.list,dims=1:20) 
ge.integrated <- IntegrateData(anchorset = ge.anchor,dims=1:20)
cat("Most Integration Complete")

saveRDS(gte.integrated, file = file.path(path_to_object, paste0(gte_filename, "_integrated.rds")) )
saveRDS(ge.integrated,  file = file.path(path_to_object, paste0(ge_filename, "_integrated.rds")) )
cat("Saved integrated seurat objects to",path_to_object,"\n")

rm("gte.list")
rm("ge.list")
rm("gte.anchors")
rm("ge.anchors")
gc()

#### =========================================== ####
#### Create Variable Feature Plots ####
#### =========================================== ####
DefaultAssay(object = gte.integrated) <- "integrated"
DefaultAssay(object = ge.integrated) <- "integrated"

# Identify the 10 most highly variable genes
size    <- 7

ifelse(!dir.exists(file.path(getwd(),parent_dir_name, "figs")),
        dir.create(file.path(getwd(),parent_dir_name, "figs"),recursive=T),
        "Directory Exists")

top10   <- head(VariableFeatures(gte.integrated), 10)
p       <- VariableFeaturePlot(gte.integrated)
p2      <- LabelPoints(plot = p, points = top10, repel = TRUE)
ggsave(file.path(path_to_object, "figs", paste0(gte_filename,"_variablegenes.tiff")), 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

top10   <- head(VariableFeatures(ge.integrated), 10)
p       <- VariableFeaturePlot(ge.integrated)
p2      <- LabelPoints(plot = p, points = top10, repel = TRUE)
ggsave(file.path(path_to_object, "figs", paste0(ge_filename,"_variablegenes.tiff")), 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

#### =========================================== ####
#### Scale, PCA, UMAP ####
#### =========================================== ####
# Standard workflow for visualization and clustering 
# (ScaleData, RunPCA, RunUMAP, FindNeighbors, FindClusters)
gte.integrated <- ScaleData(object = gte.integrated, verbose = FALSE)
ge.integrated <- ScaleData(object = ge.integrated, verbose = FALSE)

gte.integrated@meta.data$sample <- factor(gte.integrated@meta.data$sample, 
                                          levels = unique(gte.integrated@meta.data$sample))
ge.integrated@meta.data$sample <- factor(ge.integrated@meta.data$sample, 
                                          levels = unique(ge.integrated@meta.data$sample))

# Run PCA and UMAP ~20min
gte.integrated.pca <- RunPCA(object = gte.integrated, npcs = 20, verbose = FALSE)
ge.integrated.pca  <- RunPCA(object = ge.integrated, npcs = 20, verbose = FALSE)

rm("gte.integrated")
rm("ge.integrated")
gc()

gte.integrated.pca <- RunUMAP(object = gte.integrated.pca, 
                              reduction = "pca", 
                              dims = 1:20, 
                              umap.method = "uwot", 
                              metric = "cosine")

ge.integrated.pca <- RunUMAP(object = ge.integrated.pca, 
                              reduction = "pca", 
                              dims = 1:20, 
                              umap.method = "uwot", 
                              metric = "cosine")

saveRDS(gte.integrated.pca, file = file.path(path_to_object, paste0(gte_filename, "_integrated.umap.rds")) )
saveRDS(ge.integrated.pca,  file = file.path(path_to_object, paste0(ge_filename, "_integrated.umap.rds")) )
cat(paste("ScaleData, RunPCA, and RunUMAP complete.\n Saved Seurat objects to",path_to_object, "\n"))

# Determine the K-nearest neighbor graph (with first 20 PCAs)
gte.integrated.pca <- FindNeighbors(gte.integrated.pca,dims=1:20,reduction="pca")
gte.integrated.pca <- FindClusters(gte.integrated.pca, resolution = 0.3)
gte.integrated.pca <- FindClusters(gte.integrated.pca, resolution = 0.4)
saveRDS(gte.integrated.pca, file = file.path(path_to_object, paste0(gte_filename, "_integrated.umap.clustered.rds")))
cat("KNN clustering complete. Saved seurat objects to", path_to_object,"\n")

ge.integrated.pca <- FindNeighbors(ge.integrated.pca,dims=1:20,reduction="pca")
ge.integrated.pca <- FindClusters(ge.integrated.pca, resolution = 0.3)
ge.integrated.pca <- FindClusters(ge.integrated.pca, resolution = 0.4)
saveRDS(ge.integrated.pca, file = file.path(path_to_object, paste0(ge_filename, "_integrated.umap.clustered.rds")))
cat("KNN clustering complete. Saved seurat objects to", path_to_object,"\n")

#### =========================================== ####
#### Output UMAP plots by sample ####
#### =========================================== ####
size <- 5
p <- DimPlot(gte.integrated.pca, reduction = "umap", group.by = "sample") + 
    ggtitle("GRCh38+TE by Sample Origin")
ggsave(file.path(path_to_object, "figs", paste0(gte_filename,"_UMAP-sample.tiff")),
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

p <- DimPlot(ge.integrated.pca, reduction = "umap", group.by = "sample") + 
    ggtitle("GRCh38 by Sample Origin")
ggsave(file.path(path_to_object, "figs", paste0(ge_filename,"_UMAP-sample.tiff")), 
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

#### =========================================== ####
#### Output UMAP plots by cluster ####
#### =========================================== ####
p <- DimPlot(gte.integrated.pca, reduction = "umap", group.by = "cluster") + 
    ggtitle("GRCh38+TE by Sample Origin")
ggsave(file.path(path_to_object, "figs", paste0(gte_filename,"_UMAP-cluster.tiff")),
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

p <- DimPlot(ge.integrated.pca, reduction = "umap", group.by = "cluster") + 
    ggtitle("GRCh38 by Sample Origin")
ggsave(file.path(path_to_object, "figs", paste0(ge_filename,"_UMAP-cluster.tiff")),
       plot = p, units="in", width=size*3, height=size*3, dpi=300, compression = 'lzw')

cat("Seurat object processing complete. Saved to ",path_to_object,"\n")

#### End of Script ####
sessionInfo()

# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

# Matrix products: default
# BLAS/LAPACK: /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/imkl/2020.1.217/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so

# locale:
#  [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
#  [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
#  [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
#  [1] grid      stats4    parallel  stats     graphics  grDevices utils    
#  [8] datasets  methods   base     

# other attached packages:
#  [1] lawstat_3.4             ggforce_0.3.3           RColorBrewer_1.1-2     
#  [4] GO.db_3.12.1            org.Hs.eg.db_3.12.0     GOstats_2.56.0         
#  [7] graph_1.68.0            Category_2.56.0         gridExtra_2.3          
# [10] ensembldb_2.14.0        AnnotationFilter_1.14.0 GenomicFeatures_1.42.2 
# [13] AnnotationDbi_1.52.0    Biobase_2.50.0          GenomicRanges_1.42.0   
# [16] GenomeInfoDb_1.26.7     IRanges_2.24.1          S4Vectors_0.28.1       
# [19] AnnotationHub_2.22.0    BiocFileCache_1.14.0    dbplyr_2.1.0           
# [22] BiocGenerics_0.36.1     scales_1.1.1            RCurl_1.98-1.3         
# [25] forcats_0.5.1           stringr_1.4.0           dplyr_1.0.7            
# [28] purrr_0.3.4             readr_1.4.0             tidyr_1.1.3            
# [31] tibble_3.1.3            tidyverse_1.3.0         genefilter_1.72.1      
# [34] cowplot_1.1.1           ggplot2_3.3.5           Matrix_1.3-2           
# [37] SeuratObject_4.0.0      Seurat_4.0.1           

# loaded via a namespace (and not attached):
#   [1] utf8_1.2.1                    reticulate_1.20              
#   [3] tidyselect_1.1.1              RSQLite_2.2.6                
#   [5] htmlwidgets_1.5.3             BiocParallel_1.24.1          
#   [7] Rtsne_0.15                    munsell_0.5.0                
#   [9] codetools_0.2-18              ica_1.0-2                    
#  [11] pbdZMQ_0.3-5                  future_1.21.0                
#  [13] miniUI_0.1.1.1                withr_2.4.2                  
#  [15] colorspace_2.0-2              uuid_0.1-4                   
#  [17] rstudioapi_0.13               ROCR_1.0-11                  
#  [19] tensor_1.5                    listenv_0.8.0                
#  [21] Rdpack_2.1.2                  MatrixGenerics_1.2.1         
#  [23] repr_1.1.3.9000               GenomeInfoDbData_1.2.4       
#  [25] polyclip_1.10-0               farver_2.1.0                 
#  [27] bit64_4.0.5                   parallelly_1.27.0            
#  [29] vctrs_0.3.8                   generics_0.1.0               
#  [31] R6_2.5.0                      Kendall_2.2                  
#  [33] DelayedArray_0.16.3           bitops_1.0-7                 
#  [35] spatstat.utils_2.2-0          cachem_1.0.5                 
#  [37] assertthat_0.2.1              promises_1.2.0.1             
#  [39] gtable_0.3.0                  globals_0.14.0               
#  [41] goftest_1.2-2                 rlang_0.4.10                 
#  [43] splines_4.0.2                 rtracklayer_1.50.0           
#  [45] lazyeval_0.2.2                spatstat.geom_2.2-2          
#  [47] broom_0.7.5                   BiocManager_1.30.12          
#  [49] yaml_2.2.1                    reshape2_1.4.4               
#  [51] abind_1.4-5                   modelr_0.1.8                 
#  [53] backports_1.2.1               httpuv_1.6.1                 
#  [55] RBGL_1.66.0                   tools_4.0.2                  
#  [57] ellipsis_0.3.2                spatstat.core_2.3-0          
#  [59] ggridges_0.5.3                Rcpp_1.0.7                   
#  [61] plyr_1.8.6                    base64enc_0.1-3              
#  [63] progress_1.2.2                zlibbioc_1.36.0              
#  [65] prettyunits_1.1.1             rpart_4.1-15                 
#  [67] openssl_1.4.4                 deldir_0.2-10                
#  [69] pbapply_1.4-3                 zoo_1.8-9                    
#  [71] SummarizedExperiment_1.20.0   haven_2.3.1                  
#  [73] ggrepel_0.9.1                 cluster_2.1.1                
#  [75] fs_1.5.0                      magrittr_2.0.1               
#  [77] data.table_1.14.0             scattermore_0.7              
#  [79] lmtest_0.9-38                 reprex_1.0.0                 
#  [81] RANN_2.6.1                    mvtnorm_1.1-2                
#  [83] ProtGenerics_1.22.0           fitdistrplus_1.1-5           
#  [85] matrixStats_0.58.0            hms_1.0.0                    
#  [87] patchwork_1.1.1               mime_0.11                    
#  [89] evaluate_0.14                 xtable_1.8-4                 
#  [91] XML_3.99-0.6                  readxl_1.3.1                 
#  [93] compiler_4.0.2                biomaRt_2.46.3               
#  [95] KernSmooth_2.23-18            crayon_1.4.1                 
#  [97] htmltools_0.5.1.1             mgcv_1.8-34                  
#  [99] later_1.2.0                   lubridate_1.7.10             
# [101] DBI_1.1.1                     tweenr_1.0.2                 
# [103] MASS_7.3-53.1                 rappdirs_0.3.3               
# [105] boot_1.3-25                   cli_3.0.1                    
# [107] rbibutils_2.2.3               igraph_1.2.6                 
# [109] pkgconfig_2.0.3               GenomicAlignments_1.26.0     
# [111] IRdisplay_1.0.0.9000          plotly_4.9.4.1               
# [113] spatstat.sparse_2.0-0         xml2_1.3.2                   
# [115] annotate_1.68.0               XVector_0.30.0               
# [117] AnnotationForge_1.32.0        rvest_1.0.0                  
# [119] digest_0.6.27                 sctransform_0.3.2            
# [121] RcppAnnoy_0.0.19              spatstat.data_2.1-0          
# [123] Biostrings_2.58.0             cellranger_1.1.0             
# [125] leiden_0.3.9                  uwot_0.1.10                  
# [127] GSEABase_1.52.1               curl_4.3.2                   
# [129] shiny_1.6.0                   Rsamtools_2.6.0              
# [131] lifecycle_1.0.0               nlme_3.1-152                 
# [133] jsonlite_1.7.2                viridisLite_0.4.0            
# [135] askpass_1.1                   fansi_0.4.2                  
# [137] pillar_1.6.1                  lattice_0.20-41              
# [139] fastmap_1.1.0                 httr_1.4.2                   
# [141] survival_3.2-10               interactiveDisplayBase_1.28.0
# [143] glue_1.4.2                    png_0.1-7                    
# [145] Rgraphviz_2.34.0              BiocVersion_3.12.0           
# [147] bit_4.0.4                     stringi_1.7.3                
# [149] blob_1.2.1                    memoise_2.0.0                
# [151] IRkernel_1.1.1.9000           irlba_2.3.3                  
# [153] future.apply_1.7.0  