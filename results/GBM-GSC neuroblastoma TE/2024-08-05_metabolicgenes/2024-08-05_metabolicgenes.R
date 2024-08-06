#!usr/bin/env Rscript

# R version 4.2
set.seed(34)

## Import 
library(ggplot2)
library(SeuratObject)
library(Seurat)
library(dplyr)
library(gridExtra)
library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr
conflict_prefer("filter", "dplyr") ## required in %>% dplyr

set.seed(100)
print(getwd())


### Load Data
gbm <- readRDS("D:\\backup files\\2021-11-11\\Cluster\\scratch\\gete-gbm\\results\\GBM-GSC neuroblastoma TE\\2021-09-02\\ge_gbmscIntUmap-subtypes.rds")
DefaultAssay(gbm) <- "integrated"

colnames(gbm@meta.data)
head(gbm$integrated_snn_res.0.3)
head(Idents(gbm))
#setwd("D:/backup files/getegbm/GBM-TE-scRNAseq/results/GBM-GSC neuroblastoma TE/2024-08-05_metabolicgenes")

## DEG of resolution 0.3 clusters
Idents(gbm) <- "integrated_snn_res.0.3"

## genes
genes <- data.frame(gene.symbol = c("CARNS1", "CNDP1", "CNDP2", "PEPA", "NFE2L2", "NRF2", "NRF-2"))
    # CNDP2 = PEPA
    # NFE2L2 = NRF2, NRF-2
# find available genes 
present <- genes %>% filter(gene.symbol %in% rownames(gbm))
print(present)
# gene.symbol
# 1      CARNS1
# 2       CNDP1
print(match(present$gene.symbol, rownames(gbm)))
present$expr.level <- match(present$gene.symbol, rownames(gbm))
write.csv(present, "metabolic-genes_present-in-GBM-GE.csv", row.names = FALSE)

# create feature plots

size <- 5

p <- FeaturePlot(gbm, features=present$gene.symbol,
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=2
)
ggsave("figs_ge_FeaturePlot_metabolic-genes.tiff", plot = p, units="in", width=size*2.5, height=size*1, dpi=300, compression = 'lzw')

# create dot plot
p <- DotPlot(object = gbm, features=present$gene.symbol) +
    labs( x = "", y ="") +
    theme(axis.text.x = element_text(angle = 90))

ggsave("figs_ge_DotPlot_metabolic-genes.tiff", plot = p, units="in", width=size*1, height=size*1, dpi=300, compression = 'lzw')

# create dot plot with gbm subtype annotations
gbm@meta.data$gbm_subtype <- droplevels(gbm@meta.data$gbm_subtype)
gbm@meta.data$gbm_subtype <- factor(gbm@meta.data$gbm_subtype, 
                                    levels = 
                                        c("Other", 
                                          "Mesenchymal-Proneural", 
                                          "Classical-Mesenchymal", 
                                          "Classical-Proneural", 
                                          "Mesenchymal",
                                          "Classical",
                                          "Proneural"), 
                                    ordered = TRUE)
Idents(gbm) <- "gbm_subtype"

p <- DotPlot(object = gbm, features=present$gene.symbol) +
    labs( x = "", y ="") +
    theme(axis.text.x = element_text(angle = 90))

ggsave("figs_ge_DotPlot_metabolic-genes_byGBMsubtype.tiff", plot = p, units="in", width=size*1.2, height=size*0.8, dpi=300, compression = 'lzw')
ggsave("figs_ge_DotPlot_metabolic-genes_byGBMsubtype_legend.tiff", plot = p, units="in", width=size*1.2, height=size*1, dpi=300, compression = 'lzw')


# create dot plot with cell type annotations
gbm.celltypes <- read.csv("D:\\backup files\\getegbm\\GBM-TE-scRNAseq\\results\\GBM-GSC neuroblastoma TE\\2022-11-09_celltypes_gbmsc\\gbmsc_meta_celltypes.csv")
gbm@meta.data$int03_celltypes <- gbm.celltypes$int03_celltypes
rm(gbm.celltypes)
gc()
gbm@meta.data$int03_celltypes <- factor(gbm@meta.data$int03_celltypes, 
                                        levels = 
                                            c("Intermediate", 
                                              "Endothelial", 
                                              "T-cell", 
                                              "Proliferating", 
                                              "OPC",
                                              "Diff. Oligodendrocytes",
                                              "Microglia",
                                              "Astrocyte"), 
                                        ordered = TRUE)

Idents(gbm) <- "int03_celltypes"

p <- DotPlot(object = gbm, features=present$gene.symbol) +
    labs( x = "", y ="") +
    theme(axis.text.x = element_text(angle = 90))

ggsave("figs_ge_DotPlot_metabolic-genes_byCellTypes.tiff", plot = p, units="in", width=size*1.2, height=size*0.8, dpi=300, compression = 'lzw')

# create violin plots by patient
Idents(gbm) <- "int03_celltypes"
p <- VlnPlot(gbm, features=present$gene.symbol, stack=TRUE, flip=TRUE) +
    facet_grid(rows=vars(gbm@meta.data$sampleCombined), scales="free_y")+
    labs(x="")
ggsave("figs_ge_VlnPlot_metabolic-genes_byCellTypes-bySampleCombined.tiff", plot = p, units="in", width=size*1.2, height=size*1.5, dpi=300, compression = 'lzw')

Idents(gbm) <- "integrated_snn_res.0.3"
p <- VlnPlot(gbm, features=present$gene.symbol, stack=TRUE, flip=TRUE) +
    facet_grid(rows=vars(gbm@meta.data$sampleCombined), scales="free_y")+
    labs(x="")
ggsave("figs_ge_VlnPlot_metabolic-genes_byCluster-bySampleCombined.tiff", plot = p, units="in", width=size*1.2, height=size*1.5, dpi=300, compression = 'lzw')

## Sessioninfo

sessionInfo()

# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#     [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# attached base packages:
#     [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#     [1] conflicted_1.1.0   gridExtra_2.3      dplyr_1.0.10       Seurat_4.2.0      
# [5] sp_1.5-0           SeuratObject_4.1.2 ggplot2_3.3.6     
# 
# loaded via a namespace (and not attached):
#     [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6         
# [4] ellipsis_0.3.2        ggridges_0.5.4        rstudioapi_0.14      
# [7] spatstat.data_3.0-0   farver_2.1.1          leiden_0.4.3         
# [10] listenv_0.8.0         ggrepel_0.9.1         fansi_1.0.3          
# [13] codetools_0.2-18      splines_4.2.2         cachem_1.0.6         
# [16] knitr_1.40            polyclip_1.10-4       jsonlite_1.8.3       
# [19] ica_1.0-3             cluster_2.1.4         png_0.1-7            
# [22] rgeos_0.5-9           uwot_0.1.14           shiny_1.7.3          
# [25] sctransform_0.3.5     spatstat.sparse_3.0-0 compiler_4.2.2       
# [28] httr_1.4.4            assertthat_0.2.1      Matrix_1.5-1         
# [31] fastmap_1.1.0         lazyeval_0.2.2        cli_3.4.1            
# [34] later_1.3.0           htmltools_0.5.3       tools_4.2.2          
# [37] igraph_1.3.5          gtable_0.3.1          glue_1.6.2           
# [40] RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.9           
# [43] scattermore_0.8       vctrs_0.5.0           nlme_3.1-160         
# [46] progressr_0.11.0      lmtest_0.9-40         spatstat.random_3.0-1
# [49] xfun_0.34             stringr_1.4.1         globals_0.16.1       
# [52] mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.3      
# [55] irlba_2.3.5.1         goftest_1.2-3         future_1.28.0        
# [58] MASS_7.3-58.1         zoo_1.8-11            scales_1.2.1         
# [61] spatstat.core_2.4-4   promises_1.2.0.1      spatstat.utils_3.0-1 
# [64] parallel_4.2.2        RColorBrewer_1.1-3    yaml_2.3.6           
# [67] memoise_2.0.1         reticulate_1.26       pbapply_1.5-0        
# [70] rpart_4.1.19          stringi_1.7.8         rlang_1.0.6          
# [73] pkgconfig_2.0.3       matrixStats_0.62.0    evaluate_0.17        
# [76] lattice_0.20-45       ROCR_1.0-11           purrr_0.3.5          
# [79] tensor_1.5            labeling_0.4.2        patchwork_1.1.2      
# [82] htmlwidgets_1.5.4     cowplot_1.1.1         tidyselect_1.2.0     
# [85] parallelly_1.32.1     RcppAnnoy_0.0.20      plyr_1.8.7           
# [88] magrittr_2.0.3        R6_2.5.1              generics_0.1.3       
# [91] DBI_1.2.3             mgcv_1.8-41           pillar_1.8.1         
# [94] withr_2.5.0           fitdistrplus_1.1-8    survival_3.4-0       
# [97] abind_1.4-5           tibble_3.1.8          future.apply_1.9.1   
# [100] crayon_1.5.2          KernSmooth_2.23-20    utf8_1.2.2           
# [103] spatstat.geom_3.0-3   plotly_4.10.1         rmarkdown_2.17       
# [106] grid_4.2.2            data.table_1.14.4     digest_0.6.30        
# [109] xtable_1.8-4          tidyr_1.2.1           httpuv_1.6.6         
# [112] munsell_0.5.0         viridisLite_0.4.1
