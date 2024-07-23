#!usr/bin/env Rscript

# R version 4.2

## Import arguements
library(ggplot2)
library(SeuratObject)
library(Seurat)
library(dplyr)
library(gridExtra)
library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr

set.seed(100)
print(getwd())

### Load Data
gte <- readRDS("D:\\backup files\\2021-11-11\\Cluster\\scratch\\gete-gbm\\results\\GBM-GSC neuroblastoma TE\\2021-09-02\\gte_gbmscIntUmap-subtypes.rds")
DefaultAssay(gte) <- "integrated"
size = 5

colnames(gte@meta.data)
head(gte$integrated_snn_res.0.3)
head(Idents(gte))
#setwd("D:\\backup files\\getegbm\\gete-gbm\\results\\GBM-GSC neuroblastoma TE\\2024-06-25_grpathways")

## DEG of resolution 0.3 clusters
Idents(gte) <- "integrated_snn_res.0.3"

## Create Featureplots
## Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3037419/

size <- 3
gr <- c('GRIN1', 'GRIN2A', 'GRIN2B', 'GRIN3', 'GluR2'
        'GRM1', 'GRM2', 'GRM3', 'GRM4', 'GRM5', 'GRM8',
        #'mGluR1', 'mGluR2', 'mGluR3', 'mGluR4', 'mGluR5', 'mGluR8'
        'GLUD1', 'GAD1', 'GRIA1', 'GRIK2', 'SLC1A2', 'GAD1', 
        'GRIA2', 'SLC1A2', 'SLC1A3'
        )
# remove unavailable features
gr.present <- gr[-which(is.na(match(gr, rownames(gte))))] 
print(match(gr.present, rownames(gte)))

p <- FeaturePlot(gte, features=gr.present,
                 pt.size = 1, 
                 cols=c("lightyellow", "darkturquoise") 
                 )
ggsave("fig_gte_FeaturePlots_glutamate-receptors.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

## Featureplots of GENECARDS: "glutamate receptor": 1781 genes
## excluding pseudogenes and Uncategorized entries: 1767 genes 
# genecard.gr <- read.csv("GENECARD_glutamate-receptor_SEARCHRESULTS.csv")
genecard.gr <- read.csv("GENECARD_glutamate-receptor_noPseudogenes_noUncatgeorized_SEARCHRESULTS.csv")

head(genecard.gr$Gene.Symbol)
gr.present <- sort(genecard.gr$Gene.Symbol[ -which(is.na(match(genecard.gr$Gene.Symbol, rownames(gte)))) ] )
#271 genes expressed.
print(gr.present)

size <- 15

p <- FeaturePlot(gte, features=gr.present[1:50],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("fig_gte_FeaturePlots_glutamate-receptors_genecard_1-50.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

p <- FeaturePlot(gte, features=gr.present[51:100],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("fig_gte_FeaturePlots_glutamate-receptors_genecard_51-100.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

p <- FeaturePlot(gte, features=gr.present[101:150],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("fig_gte_FeaturePlots_glutamate-receptors_genecard_101-150.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

p <- FeaturePlot(gte, features=gr.present[151:200],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("fig_gte_FeaturePlots_glutamate-receptors_genecard_151-200.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

p <- FeaturePlot(gte, features=gr.present[201:250],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("fig_gte_FeaturePlots_glutamate-receptors_genecard_201-250.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

p <- FeaturePlot(gte, features=gr.present[251:length(gr.present)],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=7
)
ggsave("fig_gte_FeaturePlots_glutamate-receptors_genecard_251-.tiff", plot = p, units="in", width=size*2, height=size*1.2, dpi=300, compression = 'lzw')

## save gr.present to csv
write.csv(gr.present, file ="GENECARDS_gr-in-gte.csv", row.names=FALSE)

# Repeat the same procedure for GE dataset

### Load Data
gte <- readRDS("D:\\backup files\\2021-11-11\\Cluster\\scratch\\gete-gbm\\results\\GBM-GSC neuroblastoma TE\\2021-09-02\\ge_gbmscIntUmap-subtypes.rds")
DefaultAssay(gte) <- "integrated"

colnames(gte@meta.data)
head(gte$integrated_snn_res.0.3)
head(Idents(gte))
#setwd("D:\\backup files\\getegbm\\gete-gbm\\results\\GBM-GSC neuroblastoma TE\\2024-06-25_grpathways")

## DEG of resolution 0.3 clusters
Idents(gte) <- "integrated_snn_res.0.3"

## Featureplots of GENECARDS: "glutamate receptor": 1781 genes
## excluding pseudogenes and Uncategorized entries: 1767 genes 
# genecard.gr <- read.csv("GENECARD_glutamate-receptor_SEARCHRESULTS.csv")
genecard.gr <- read.csv("GENECARD_glutamate-receptor_noPseudogenes_noUncatgeorized_SEARCHRESULTS.csv")

head(genecard.gr$Gene.Symbol)
gr.present <- sort(genecard.gr$Gene.Symbol[ -which(is.na(match(genecard.gr$Gene.Symbol, rownames(gte)))) ] )
#271 genes expressed.
print(gr.present)

size <- 15

p <- FeaturePlot(gte, features=gr.present[1:50],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("figs_ge_FeaturePlots_glutamate-receptors_genecard_1-50.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

p <- FeaturePlot(gte, features=gr.present[51:100],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("figs_ge_FeaturePlots_glutamate-receptors_genecard_51-100.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

p <- FeaturePlot(gte, features=gr.present[101:150],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("figs_ge_FeaturePlots_glutamate-receptors_genecard_101-150.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

p <- FeaturePlot(gte, features=gr.present[151:200],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("figs_ge_FeaturePlots_glutamate-receptors_genecard_151-200.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

p <- FeaturePlot(gte, features=gr.present[201:250],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=8
)
ggsave("figs_ge_FeaturePlots_glutamate-receptors_genecard_201-250.tiff", plot = p, units="in", width=size*2, height=size*2, dpi=300, compression = 'lzw')

p <- FeaturePlot(gte, features=gr.present[251:length(gr.present)],
                 pt.size = 1, 
                 cols=c("lightyellow", "steelblue3"),
                 ncol=7
)
ggsave("figs_ge_FeaturePlots_glutamate-receptors_genecard_251-.tiff", plot = p, units="in", width=size*2, height=size*1.2, dpi=300, compression = 'lzw')

## save gr.present to csv
write.csv(gr.present, file ="GENECARDS_gr-in-ge.csv", row.names=FALSE)



## Session Info
sessionInfo()

# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#     [1] LC_COLLATE=English_United States.utf8 
# [2] LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# attached base packages:
#     [1] stats     graphics  grDevices utils     datasets  methods  
# [7] base     
# 
# other attached packages:
#     [1] gridExtra_2.3      conflicted_1.1.0   dplyr_1.0.10      
# [4] Seurat_4.2.0       sp_1.5-0           SeuratObject_4.1.2
# [7] ggplot2_3.3.6     
# 
# loaded via a namespace (and not attached):
#     [1] Rtsne_0.16            colorspace_2.0-3     
# [3] deldir_1.0-6          ellipsis_0.3.2       
# [5] ggridges_0.5.4        rstudioapi_0.14      
# [7] spatstat.data_3.0-0   farver_2.1.1         
# [9] leiden_0.4.3          listenv_0.8.0        
# [11] ggrepel_0.9.1         fansi_1.0.3          
# [13] codetools_0.2-18      splines_4.2.2        
# [15] cachem_1.0.6          knitr_1.40           
# [17] polyclip_1.10-4       jsonlite_1.8.3       
# [19] ica_1.0-3             cluster_2.1.4        
# [21] png_0.1-7             rgeos_0.5-9          
# [23] uwot_0.1.14           shiny_1.7.3          
# [25] sctransform_0.3.5     spatstat.sparse_3.0-0
# [27] compiler_4.2.2        httr_1.4.4           
# [29] assertthat_0.2.1      Matrix_1.5-1         
# [31] fastmap_1.1.0         lazyeval_0.2.2       
# [33] cli_3.4.1             later_1.3.0          
# [35] htmltools_0.5.3       tools_4.2.2          
# [37] igraph_1.3.5          gtable_0.3.1         
# [39] glue_1.6.2            RANN_2.6.1           
# [41] reshape2_1.4.4        Rcpp_1.0.9           
# [43] scattermore_0.8       vctrs_0.5.0          
# [45] nlme_3.1-160          progressr_0.11.0     
# [47] lmtest_0.9-40         spatstat.random_3.0-1
# [49] xfun_0.34             stringr_1.4.1        
# [51] globals_0.16.1        mime_0.12            
# [53] miniUI_0.1.1.1        lifecycle_1.0.3      
# [55] irlba_2.3.5.1         goftest_1.2-3        
# [57] future_1.28.0         MASS_7.3-58.1        
# [59] zoo_1.8-11            scales_1.2.1         
# [61] spatstat.core_2.4-4   promises_1.2.0.1     
# [63] spatstat.utils_3.0-1  parallel_4.2.2       
# [65] RColorBrewer_1.1-3    yaml_2.3.6           
# [67] memoise_2.0.1         reticulate_1.26      
# [69] pbapply_1.5-0         rpart_4.1.19         
# [71] stringi_1.7.8         rlang_1.0.6          
# [73] pkgconfig_2.0.3       matrixStats_0.62.0   
# [75] evaluate_0.17         lattice_0.20-45      
# [77] ROCR_1.0-11           purrr_0.3.5          
# [79] tensor_1.5            labeling_0.4.2       
# [81] patchwork_1.1.2       htmlwidgets_1.5.4    
# [83] cowplot_1.1.1         tidyselect_1.2.0     
# [85] parallelly_1.32.1     RcppAnnoy_0.0.20     
# [87] plyr_1.8.7            magrittr_2.0.3       
# [89] R6_2.5.1              generics_0.1.3       
# [91] DBI_1.1.3             mgcv_1.8-41          
# [93] pillar_1.8.1          withr_2.5.0          
# [95] fitdistrplus_1.1-8    survival_3.4-0       
# [97] abind_1.4-5           tibble_3.1.8         
# [99] future.apply_1.9.1    crayon_1.5.2         
# [101] KernSmooth_2.23-20    utf8_1.2.2           
# [103] spatstat.geom_3.0-3   plotly_4.10.1        
# [105] rmarkdown_2.17        grid_4.2.2           
# [107] data.table_1.14.4     digest_0.6.30        
# [109] xtable_1.8-4          tidyr_1.2.1          
# [111] httpuv_1.6.6          munsell_0.5.0        
# [113] viridisLite_0.4.1


