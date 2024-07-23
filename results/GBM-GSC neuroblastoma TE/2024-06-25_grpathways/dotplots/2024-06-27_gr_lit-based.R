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
#setwd("D:\\backup files\\getegbm\\gete-gbm\\results\\GBM-GSC neuroblastoma TE\\2024-06-25_grpathways")

## DEG of resolution 0.3 clusters
Idents(gbm) <- "integrated_snn_res.0.3"

## glutamate associated genes
gr.genes <- read.csv("glutamate-genes.csv", header = TRUE)
head(gr.genes$gene.symbol)
# find available gr.genes 
gr.present <- gr.genes %>% filter(gene.symbol %in% rownames(gbm))
print(gr.present)
print(match(gr.present$gene.symbol, rownames(gbm)))
gr.present$expr.level <- match(gr.present$gene.symbol, rownames(gbm))
write.csv(gr.present, "glutamate-genes_present-in-GBM-GE.csv", row.names = FALSE)

# create feature plots

size <- 5

p <- FeaturePlot(gbm, features=gr.present$gene.symbol,
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=3
)
ggsave("figs_ge_FeaturePlot_glutamate-receptors.tiff", plot = p, units="in", width=size*3.5, height=size*2, dpi=300, compression = 'lzw')

# create dot plot
p <- DotPlot(object = gbm, features=gr.present$gene.symbol) +
    labs( x = "", y ="") +
    theme(axis.text.x = element_text(angle = 90))

ggsave("figs_ge_DotPlot_glutamate-receptors.tiff", plot = p, units="in", width=size*1, height=size*1, dpi=300, compression = 'lzw')

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

p <- DotPlot(object = gbm, features=gr.present$gene.symbol) +
    labs( x = "", y ="") +
    theme(axis.text.x = element_text(angle = 90))

ggsave("figs_ge_DotPlot_glutamate-receptors_byGBMsubtype.tiff", plot = p, units="in", width=size*1.2, height=size*0.8, dpi=300, compression = 'lzw')


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

p <- DotPlot(object = gbm, features=gr.present$gene.symbol) +
    labs( x = "", y ="") +
    theme(axis.text.x = element_text(angle = 90))

ggsave("figs_ge_DotPlot_glutamate-receptors_byCellTypes.tiff", plot = p, units="in", width=size*1.2, height=size*0.8, dpi=300, compression = 'lzw')


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
#     [1] conflicted_1.1.0   gridExtra_2.3      dplyr_1.0.10       Seurat_4.2.0       sp_1.5-0          
# [6] SeuratObject_4.1.2 ggplot2_3.3.6     
# 
# loaded via a namespace (and not attached):
#     [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6          ellipsis_0.3.2       
# [5] ggridges_0.5.4        rstudioapi_0.14       spatstat.data_3.0-0   leiden_0.4.3         
# [9] listenv_0.8.0         farver_2.1.1          ggrepel_0.9.1         fansi_1.0.3          
# [13] codetools_0.2-18      splines_4.2.2         cachem_1.0.6          polyclip_1.10-4      
# [17] jsonlite_1.8.3        ica_1.0-3             cluster_2.1.4         png_0.1-7            
# [21] rgeos_0.5-9           uwot_0.1.14           shiny_1.7.3           sctransform_0.3.5    
# [25] spatstat.sparse_3.0-0 compiler_4.2.2        httr_1.4.4            assertthat_0.2.1     
# [29] Matrix_1.5-1          fastmap_1.1.0         lazyeval_0.2.2        cli_3.4.1            
# [33] later_1.3.0           htmltools_0.5.3       tools_4.2.2           igraph_1.3.5         
# [37] gtable_0.3.1          glue_1.6.2            RANN_2.6.1            reshape2_1.4.4       
# [41] Rcpp_1.0.9            scattermore_0.8       vctrs_0.5.0           nlme_3.1-160         
# [45] progressr_0.11.0      lmtest_0.9-40         spatstat.random_3.0-1 stringr_1.4.1        
# [49] globals_0.16.1        mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.3      
# [53] irlba_2.3.5.1         goftest_1.2-3         future_1.28.0         MASS_7.3-58.1        
# [57] zoo_1.8-11            scales_1.2.1          spatstat.core_2.4-4   promises_1.2.0.1     
# [61] spatstat.utils_3.0-1  parallel_4.2.2        RColorBrewer_1.1-3    memoise_2.0.1        
# [65] reticulate_1.26       pbapply_1.5-0         rpart_4.1.19          stringi_1.7.8        
# [69] rlang_1.0.6           pkgconfig_2.0.3       matrixStats_0.62.0    lattice_0.20-45      
# [73] ROCR_1.0-11           purrr_0.3.5           tensor_1.5            patchwork_1.1.2      
# [77] htmlwidgets_1.5.4     labeling_0.4.2        cowplot_1.1.1         tidyselect_1.2.0     
# [81] parallelly_1.32.1     RcppAnnoy_0.0.20      plyr_1.8.7            magrittr_2.0.3       
# [85] R6_2.5.1              generics_0.1.3        DBI_1.2.3             pillar_1.8.1         
# [89] withr_2.5.0           mgcv_1.8-41           fitdistrplus_1.1-8    survival_3.4-0       
# [93] abind_1.4-5           tibble_3.1.8          future.apply_1.9.1    crayon_1.5.2         
# [97] KernSmooth_2.23-20    utf8_1.2.2            spatstat.geom_3.0-3   plotly_4.10.1        
# [101] grid_4.2.2            data.table_1.14.4     digest_0.6.30         xtable_1.8-4         
# [105] tidyr_1.2.1           httpuv_1.6.6          munsell_0.5.0         viridisLite_0.4.1 
