Author: Samantha Yuen
Date: 2024-06-26 to 2024-06-25

# GeneCards -- Search for "glutamate receptor" 

- Obtained 1781 genes
- Further filtered out pseudogenes and uncategorized entries resulting in 1767 genes
- GENECARD_glutamate-receptor_noPseudogenes_noUncatgeorized_SEARCHRESULTS.csv
- Older literature to test the initial sub-20 genes earlier in the code: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3421042/ (Prickett & Samuels, 2013. Clin Cancer Res)
- Created FeaturePlots of resulting 271 genes present in scRNA-seq human GBM assay (50 plots per TIFF figure). 

# Gene Ontology -- Search results for "glutamate receptor" on AmiGO2

`AmiGO_glutamate-receptor_SEARCHRESULTS.tsv`
- Interpretation: https://www.researchgate.net/figure/Gene-Ontology-GO-analysis-of-biopsy-microarray-data-A-GO-analysis-of-biological_fig1_305268803

GE
- 3841 significantly differentially expressed genes our of 7566 variable genes
- 80137 GOIDs across 12 clusters resulted from a non-biased analysis
- 153 Terms across 12 clusters (0-11) were related to glutamate receptor terms (AmiGO)
  - (53 glutamate receptor GO Terms * 12 clusters = 636 hits maximum)
- only 1 TERM significantly enriched.
- (≥1.25-fold, P ≤ 0.05); Only [GO:0098990 AMPA selective glutamate receptor signaling pathway] seems to be significant across all clusters (mainly cluster 0 and 1), but this GOterm only has 1 gene in total, so it's not indicative of glutamate receptor expression activation. 
- 1 term was significantly enriched 
```
> nrow(deg)
[1] 7566
> deg <- dplyr::filter(deg, p_val_adj<=0.05 & 
+                   (avg_log2FC>=0.5 | avg_log2FC<=-0.5))
> nrow(deg)
[1] 3841
> glimpse(deg)
Rows: 3,841
Columns: 7

> deg.go <- read.csv("ge_GO_GR.csv", header = TRUE, row.names = 1)
> dim(deg.go)
[1] 153  11

> deg.go <- read.csv("ge_GO_GR.csv", header = TRUE, row.names = 1)
> deg.go <- deg.go %>% dplyr::filter(p.adj <= 0.05) %>% dplyr::group_by(cluster)
> deg.go
# A tibble: 1 × 11
# Groups:   cluster [1]
  GOID       Term                      Class Pvalue  p.adj nlog1…¹ OddsR…² ExpCo…³ Count  Size cluster
  <chr>      <chr>                     <chr>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl> <int> <int> <fct>  
1 GO:0098990 AMPA selective glutamate… BP    0.0115 0.0413    1.38     Inf  0.0115     1     1 1      
# … with abbreviated variable names ¹​nlog10p.adj, ²​OddsRatio, ³​ExpCount
```

GTE:

```
> nrow(deg)
[1] 5656
> deg <- dplyr::filter(deg, p_val_adj<=0.05 & 
+                          (avg_log2FC>=0.5 | avg_log2FC<=-0.5))
> nrow(deg)
[1] 2815
...
> summary(deg.go$p.adj <0.05)
   Mode   FALSE    TRUE 
logical   49113   16248 
```
```
# A tibble: 1 × 11
# Groups:   cluster [1]
  GOID       Class Pvalue  p.adj    nlog10p.adj, OddsRatio, ExpCount Count  Size    cluster
  <chr>      <chr>  <dbl>  <dbl>    <dbl>        <dbl>      <dbl>    <int> <int>    <int>
1 GO:0098990 BP    0.0115 0.0413    1.38         Inf        0.0115   1      1       1
```

# KEGG_grpathways.csv -- Search results for "glutamate receptor" 

- https://www.kegg.jp/kegg/pathway.html 
- 17 hits were obtained from KEGG database related to glutamate receptors
  - KEGG entries and their codes were saved in a csv file. 

Results GE:
- Obtained 1368 KEGG Terms with 5 clusters that are signficantly enirched in one glutmate receptor related term 
- Around <10% [7-14 genes] of genes in our assay are present in KEGGID 04540 which is better than what we can say from the GO analysis.

```
# hsa04540    Gap junction -- https://www.kegg.jp/entry/hsa04540 Cellular Process
# # A tibble: 5 × 11
# # Groups:   cluster [5]
# KEGGID Term  Class     Pvalue    p.adj nlog10p.adj OddsRatio ExpCount Count  Size cluster
# <chr>  <chr> <chr>      <dbl>    <dbl>       <dbl>     <dbl>    <dbl> <int> <int> <chr>  
#     1 04540  NA    KEGG  0.000750   0.0115          1.94      4.56     1.98     8    90 2      
# 2 04540  NA    KEGG  0.0000169  0.000193        3.71      6.21     1.90    10    90 7      
# 3 04540  NA    KEGG  0.000329   0.00738         2.13      3.90     3.22    11    90 8      
# 4 04540  NA    KEGG  0.00125    0.00740         2.13      4.74     1.66     7    90 9      
# 5 04540  NA    KEGG  0.00000544 0.000171        3.77      4.96     3.39    14    90 11 
```

GTE dataset
- 1368 KEGGIDs obtained from analysis, but 0 were significantly enriched for glutamate receptor related terms. 

# Extract GeneIDs per GO/KEGG ID (GE)

```
> genes.df <- allGOgenes(deg.go)
Processing cluster 1
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 5069
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 706
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 507
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 125
Genes Extracted for cluster: 1
Processing cluster 2
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 4981
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 686
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 447
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 115
Genes Extracted for cluster: 2
Processing cluster 3
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 5248
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 786
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 544
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 144
Genes Extracted for cluster: 3
Processing cluster 4
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 2999
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 371
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 335
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 60
Genes Extracted for cluster: 4
Processing cluster 5
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 5475
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 782
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 603
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 136
Genes Extracted for cluster: 5
Processing cluster 6
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 7176
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 1162
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 647
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 169
Genes Extracted for cluster: 6
Processing cluster 7
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 4501
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 670
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 463
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 117
Genes Extracted for cluster: 7
Processing cluster 8
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 5322
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 725
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 574
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 118
Genes Extracted for cluster: 8
Processing cluster 9
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 6729
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 1073
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 673
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 157
Genes Extracted for cluster: 9
Processing cluster 10
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 4827
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 680
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 485
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 118
Genes Extracted for cluster: 10
Processing cluster 11
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 5290
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 750
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 553
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 132
Genes Extracted for cluster: 11
Processing cluster 12
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 6574
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 1105
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 655
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 157
Genes Extracted for cluster: 12
Genes Extracted for all clusters
```

# Extract GeneIDs per GO/KEGG ID (GTE)

```
> deg.go <- readRDS("gte_GO-KEGG.rds")
> genes.df <- allGOgenes(deg.go)
Processing cluster 1
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 4478
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 597
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 456
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 121
Genes Extracted for cluster: 1
Processing cluster 2
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 4728
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 612
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 492
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 132
Genes Extracted for cluster: 2
Processing cluster 3
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 4313
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 573
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 420
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 116
Genes Extracted for cluster: 3
Processing cluster 4
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 2259
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 267
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 263
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 59
Genes Extracted for cluster: 4
Processing cluster 5
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 6740
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 1063
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 638
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 157
Genes Extracted for cluster: 5
Processing cluster 6
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 4816
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 701
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 551
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 138
Genes Extracted for cluster: 6
Processing cluster 7
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 5674
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 903
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 576
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 146
Genes Extracted for cluster: 7
Processing cluster 8
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 4097
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 586
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 499
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 108
Genes Extracted for cluster: 8
Processing cluster 9
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 5008
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 678
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 577
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 122
Genes Extracted for cluster: 9
Processing cluster 10
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 4469
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 614
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 467
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 117
Genes Extracted for cluster: 10
Processing cluster 11
Processing Ontology Category: BP
Number of GO terms in this Ontology Class: 5810
Processing Ontology Category: MF
Number of GO terms in this Ontology Class: 904
Processing Ontology Category: CC
Number of GO terms in this Ontology Class: 565
Processing Ontology Category: KEGG
Number of GO terms in this Ontology Class: 152
Genes Extracted for cluster: 11
Genes Extracted for all clusters
```

# Environment

Sourced util_go.R from /bin.

Environment: 
R 4.2.0
Seurat v4
gridExtra 2.3
ggplot2 3.36
tidyverse 1.3.2
Bioconductor 3.16
    GOstats 2.64
    KEGGREST 1.38.0 (replaces KEGG.db)
