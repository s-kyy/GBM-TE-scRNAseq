#!usr/bin/env Rscript

# R version 4.2

set.seed(34)

## Import
library(ggplot2)
library(tidyverse)
library(dplyr)
library(conflicted)
conflict_prefer("select", "dplyr") ## required in %>% dplyr

set.seed(100)
print(getwd())

source('../../../../bin/util_go.R') 

## Gene Ontology and KEGG analysis of GE
writeLines("Loading differentially expressed genes csv")
deg <- read.csv("D:\\backup files\\getegbm\\GBM-TE-scRNAseq\\results\\GBM-GSC neuroblastoma TE\\2023-07-12_DEG\\2023-08-08\\ge_markers_0.3.csv", row.names = 1)
rownames(deg) <- NULL
# deg <- deg %>% tibble::rownames_to_column(var = "gene")

# Ensure cluster column contains integers (not factors)
if (sapply(deg,mode)[length(colnames(deg))-1] == "numeric'") {
    writeLines("Cluster column contains integers")
    continue
} else {
    writeLines("Converting cluster column to integers")
    deg$cluster <- as.integer(deg$cluster)
}

# Filter non-significant genes
nrow(deg)
deg <- dplyr::filter(deg, p_val_adj<=0.05 & 
                  (avg_log2FC>=0.5 | avg_log2FC<=-0.5))
nrow(deg)
glimpse(deg)

deg.go <- geneontology(deg)
saveRDS(deg.go, file="ge_GO-KEGG.rds")

deg.go <- makeGOsummary(deg.go)
write.csv(deg.go, "ge_GO_summary.csv")

head(deg.go)

summary(deg.go$p.adj <0.05)
savehistory("2024-06-26_ge_GO.Rhistory")

### Analyze GOID for Glutamate Receptors
# deg.go <- read.csv("ge_GO_summary.csv", row.names = 1)
# rownames(deg.go) <- NULL
amigo.gr <- read.delim("AmiGO_glutamate-receptor_SEARCHRESULTS.tsv", header = FALSE)
colnames(amigo.gr) <- c(colnames(deg.go)[1:3], "Synonyms")

glimpse(deg.go)

deg.go <- deg.go %>% dplyr::semi_join(amigo.gr, by = "GOID")
glimpse(deg.go)
write.csv(deg.go, "ge_GO_GR.csv")
deg.go <- read.csv("ge_GO_GR.csv", header = TRUE, row.names = 1)

deg.go <- deg.go %>% dplyr::filter(p.adj <= 0.05) %>% dplyr::group_by(cluster)

size <- 10
#deg.go$cluster <- factor(deg.go$cluster)
# p <- ggplot(deg.go, aes(x=nlog10p.adj, y=Term, fill=Class)) + 
#     geom_bar(stat='identity') + xlab("-log10(p)") + ylab("GO Term") +
#     scale_fill_manual("GO Class", values = c("BP" = "slateblue3", "MF" = "lightseagreen", "CC" = "darkgoldenrod1")) +
#     geom_text(aes(label=Count), color="azure4", hjust = -0.5) + theme_classic()
#     # facet_wrap(~cluster)
# ggsave("fig_ge_GO_barplot_glutamate-receptors.tiff", plot = p, units="in", width=size*5, height=size*2, dpi=300, compression = 'lzw')

### Analyze KEGGIDs for Glutamate Receptors
kegg.gr <- read.csv("KEGG_grpathways.csv", header = TRUE)
kegg.gr <- kegg.gr %>% mutate(KEGGID= str_extract(KEGG_entry, "[0-9]{5}$") %>% 
                       as.character())
deg.go <- readRDS("ge_GO-KEGG.rds")
deg.kegg <- makeKEGGsummary(deg.go)
write.csv(deg.kegg, "ge_KEGG_summary.csv")

deg.kegg <- deg.kegg %>% dplyr::semi_join(kegg.gr, by = "KEGGID")
glimpse(deg.kegg)
write.csv(deg.kegg, "ge_KEGG_GR.csv")
# deg.kegg <- read.csv("ge_KEGG_GR.csv", header = TRUE, row.names = 1)

deg.kegg %>% dplyr::filter(p.adj <=0.05) %>% group_by(cluster)
# # A tibble: 5 × 11
# # Groups:   cluster [5]
# KEGGID Term  Class     Pvalue    p.adj nlog10p.adj OddsRatio ExpCount Count  Size cluster
# <chr>  <chr> <chr>      <dbl>    <dbl>       <dbl>     <dbl>    <dbl> <int> <int> <chr>  
#     1 04540  NA    KEGG  0.000750   0.0115          1.94      4.56     1.98     8    90 2      
# 2 04540  NA    KEGG  0.0000169  0.000193        3.71      6.21     1.90    10    90 7      
# 3 04540  NA    KEGG  0.000329   0.00738         2.13      3.90     3.22    11    90 8      
# 4 04540  NA    KEGG  0.00125    0.00740         2.13      4.74     1.66     7    90 9      
# 5 04540  NA    KEGG  0.00000544 0.000171        3.77      4.96     3.39    14    90 11 

gc()

### Export genesbyID per cluster

deg.go <- readRDS("ge_GO-KEGG.rds")
# deg.go.c1 <- deg.go$`0`[[no.classes+1]]
# rm(deg.go)
# gc()
# summary(deg.go.c1)
# genes.list <- genesbyGO(deg.go.c1,1)

genes.df <- allGOgenes(deg.go)
write.csv(genes.df, "ge_geneIDbyGO-KEGG.csv", row.names = FALSE)
gc()

colnames(genes.df)
head(genes.df)
genes.df <- as.data.frame(genes.df)
genes.df %>% dplyr::filter(ID == "GO:0098990" & cluster == "1")
# ID Genes Class cluster
# GO.0098990.1 GO:0098990  PLP1    BP       1

### Repeat for GTE analysis. 
# restart R
writeLines("Loading differentially expressed genes csv")
deg <- read.csv("D:\\backup files\\getegbm\\GBM-TE-scRNAseq\\results\\GBM-GSC neuroblastoma TE\\2023-07-12_DEG\\2023-08-08\\gte_markers_0.3.csv", row.names = 1)
rownames(deg) <- NULL
# deg <- deg %>% tibble::rownames_to_column(var = "gene")

# Ensure cluster column contains integers (not factors)
if (sapply(deg,mode)[length(colnames(deg))-1] == "numeric'") {
    writeLines("Cluster column contains integers")
    continue
} else {
    writeLines("Converting cluster column to integers")
    deg$cluster <- as.integer(deg$cluster)
}

# Filter non-significant genes
nrow(deg)
deg <- dplyr::filter(deg, p_val_adj<=0.05 & 
                         (avg_log2FC>=0.5 | avg_log2FC<=-0.5))
nrow(deg)
glimpse(deg)

## Perform GO + KEGG Analysis
deg.go <- geneontology(deg)
saveRDS(deg.go, file="gte_GO-KEGG.rds")
rm(deg)
gc()

deg.go <- makeGOsummary(deg.go)
write.csv(deg.go, "gte_GO_summary.csv")

head(deg.go)

summary(deg.go$p.adj <0.05)
savehistory("2024-06-26_gte_GO.Rhistory")

### Analyze GOID for Glutamate Receptors
amigo.gr <- read.delim("AmiGO_glutamate-receptor_SEARCHRESULTS.tsv", header = FALSE)
colnames(amigo.gr) <- c(colnames(deg.go)[1:3], "Synonyms")

glimpse(deg.go)

deg.go <- deg.go %>% dplyr::semi_join(amigo.gr, by = "GOID")
glimpse(deg.go)
write.csv(deg.go, "gte_GO_GR.csv")
# deg.go <- read.csv("ge_GO_GR.csv", header = TRUE, row.names = 1)

deg.go <- deg.go %>% dplyr::filter(p.adj <= 0.05) %>% dplyr::group_by(cluster)
deg.go$cluster <- factor(deg.go$cluster)

size <- 10
# p <- ggplot(deg.go, aes(x=nlog10p.adj, y=Term, fill=Class)) + 
#     geom_bar(stat='identity') + xlab("-log10(p)") + ylab("GO Term") +
#     scale_fill_manual("GO Class", values = c("BP" = "slateblue3", "MF" = "lightseagreen", "CC" = "darkgoldenrod1")) +
#     geom_text(aes(label=Count), color="azure4", hjust = -0.5) + theme_classic()
#     # facet_wrap(~cluster)
# ggsave("fig_ge_GO_barplot_glutamate-receptors.tiff", plot = p, units="in", width=size*5, height=size*2, dpi=300, compression = 'lzw')

### Analyze KEGGIDs for Glutamate Receptors
kegg.gr <- read.csv("KEGG_grpathways.csv", header = TRUE)
kegg.gr <- kegg.gr %>% mutate(KEGGID= str_extract(KEGG_entry, "[0-9]{5}$") %>% 
                                  as.character())
deg.go <- readRDS("gte_GO-KEGG.rds")
deg.kegg <- makeKEGGsummary(deg.go)
write.csv(deg.kegg, "gte_KEGG_summary.csv")
deg.kegg <- read.csv("gte_KEGG_summary.csv", row.names = 1)
rownames(deg.kegg) <- NULL
deg.kegg$KEGGID <- as.character(deg.kegg$KEGGID)

deg.kegg <- deg.kegg %>% dplyr::semi_join(kegg.gr, by = "KEGGID")
glimpse(deg.kegg)
write.csv(deg.kegg, "gte_KEGG_GR.csv")
# deg.kegg <- read.csv("gte_KEGG_GR.csv", header = TRUE, row.names = 1)

deg.kegg %>% dplyr::filter(p.adj <=0.05) %>% group_by(cluster)

### Obtain list of GO genes

deg.go <- readRDS("gte_GO-KEGG.rds")
genes.df <- allGOgenes(deg.go)
write.csv(genes.df, "gte_geneIDbyGO-KEGG.csv", row.names = FALSE)
gc()

### Session Info
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
#     [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#     [1] GO.db_3.16.0         org.Hs.eg.db_3.16.0  KEGGREST_1.38.0      GOstats_2.64.0       graph_1.76.0        
# [6] Category_2.64.0      Matrix_1.5-1         AnnotationDbi_1.60.2 IRanges_2.32.0       S4Vectors_0.36.2    
# [11] Biobase_2.58.0       BiocGenerics_0.44.0  conflicted_1.1.0     dplyr_1.0.10         ggplot2_3.3.6       
# 
# loaded via a namespace (and not attached):
#     [1] Rcpp_1.0.9             lattice_0.20-45        png_0.1-7              Biostrings_2.66.0     
# [5] digest_0.6.30          assertthat_0.2.1       utf8_1.2.2             R6_2.5.1              
# [9] GenomeInfoDb_1.34.9    RSQLite_2.3.6          httr_1.4.4             pillar_1.8.1          
# [13] zlibbioc_1.44.0        rlang_1.0.6            rstudioapi_0.14        annotate_1.76.0       
# [17] Rgraphviz_2.42.0       blob_1.2.3             labeling_0.4.2         splines_4.2.2         
# [21] RCurl_1.98-1.14        bit_4.0.4              munsell_0.5.0          compiler_4.2.2        
# [25] pkgconfig_2.0.3        tidyselect_1.2.0       tibble_3.1.8           GenomeInfoDbData_1.2.9
# [29] XML_3.99-0.16.1        AnnotationForge_1.40.2 fansi_1.0.3            crayon_1.5.2          
# [33] withr_2.5.0            bitops_1.0-7           grid_4.2.2             RBGL_1.74.0           
# [37] xtable_1.8-4           GSEABase_1.60.0        gtable_0.3.1           lifecycle_1.0.3       
# [41] DBI_1.2.3              magrittr_2.0.3         scales_1.2.1           cli_3.4.1             
# [45] cachem_1.0.6           farver_2.1.1           XVector_0.38.0         genefilter_1.80.3     
# [49] generics_0.1.3         vctrs_0.5.0            tools_4.2.2            bit64_4.0.5           
# [53] glue_1.6.2             fastmap_1.1.0          survival_3.4-0         colorspace_2.0-3      
# [57] memoise_2.0.1
