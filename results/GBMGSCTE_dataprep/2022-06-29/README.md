## Folder contents

### Badhuri et al. 2020 GBM dataset (archive): 

Dir: GBM_rxiv
Script: seurat_brain-dataexpo.ipynb
Associated results in backup:

- 2021-03-30: QC figures before cell filtering
- 2021-03-31: QC figures after cell filtering
- 2021-04-01: UMAP w resolution 0.8
- 2021-04-06: UMAP with GBM subtype markers, QC metrics and Sample Distribution
- 2021-04-07: DEG .RData files & Exploring cell markers 
- ~~2021-04-19~~: GO results objects (by cluster and all cells combined)
- 2021-04-22 to 23: Redid the GO analysis
- 2021-04-26 to 27: Producing the GO tables (refer to 27th)
- GBM_GSEA_may21: GSEA results

## GBM + Wang et al. 2020 Neuroblastoma (GSC) Dataset:

Dir: GBMGSCTE and GBMGSC_monocle
Script: neuroblastoma.ipynb and neurolastoma2.ipynb
Associated results in backup:

- 2021-06-10: Unfiltered Seurat Object for neuroblastoma dataset
- 2021-06-11: Filtered neuroblastoma dataset + merged with GBM dataset. ge_gbmsc and gte_gbmsc. Quality control (before filtering)
- 2021-06-15: PCA and PC figures
- 2021-06-17: ~~UMAP 0.7 resolution (with updated sample names)~~, Exploring cellular markers, and matching cluster labels
- 2021-06-21: DEG results table, Cluster Analysis, Gene Ontology, table of GO results & conversion back to gene symbols, GO Barplots
- 2021-07-02: DEG and GO of GE: 19/17, and GETE: 18/14 clusters to identify the clusters, GO Barplots
- 2021-07-06: DEG and GO of GTE: 3, 10, 19 to identify marker genes, GO Barplots
- 2021-07-08: Added MGMT metadata to ge and gte Seurat objects (UMAP 0.7), UMAP of MGMT+ verus MGMT- cellS. MGMT, GLAST(SLC1A3) & CD133 (PROM1) expression grouped by clusters and sample_IDs visualized with pie charts and bar plots. 
- 2021-07-16: 
    - UMAPS showing expression levels of 'VEGFA', 'TET1', 'ESRRB', 'XRN2', 'SOD1' cancer stem cell master regulators.
    - Pie Chart of CD133, GLAST and MGMT based on cluster composition (with cluster labels).
    - Barplots and Pie Charts of CSC master regulator expression based on cluster composition and sample composition. 
    - Barplots and Pie Charts of Genes unique in human organoids based on cluster composition and sample composition.
- 2021-07-22: Percentage of TE UMI & Genes in each cluster 
- 2021-07-23: Trajectory Analysis data + preliminary figure for biweekly meeting (2021-07-23).
- GBMSC_GSEA-jul09: GSEA Results (Bar plots, UMAPs)
    - 2021-09-03: Create Barplot with percentage of GBM subtype in each sample. 
- 2021-07-27: 
    - Updated seurat objects to have GBM subtype classifications based on Verhaak and Wang gene sets
    - Variable Feature plots by sample for ge and gte. 
    - New integrated seurat objects for **neuroblastoma2.ipynb** (gbmscGEint.rds & gbmscGTEint.rds)\
    - PC evaluation figures: Heatmaps, Elbow Plot, JackStrawPlots
    - QC: UMI, Feature, & mitRatio ViolinPlots & Scatter plots
    - (To Do) Redid UMAPs 
- 2021-08-02:
    - QC: UMI, Feature & mitoRatio prior to integration
    - gbmsc_ge and gbmsc_gte: Version 2 seurat objects with SF11215 sample removed and more agressive filtering of features when creating the seurat object (min features = 1000, min cells per gene = 30).
    - QC_v2: Redid for the new seurat object
    - Find Varaible Feature plots by sample for ge and gte
    - Re-integrate seurat objects in neuroblastoma 2.
    - Plots to evaluate PC analysis (2021-08-03)
    - UMAP & Cluster Analysis (visualized a range of PCs and resolutions in UMAP space) (2021-08-03)
    - Trajectory Analysis (2021-08-03 - 2021-08-04)
        - Seurat Clusters: 7 (GE & GTE)
        - Monocle Clusters: 6 (GE) & 7 (GTE)
- 2021-08-04: Performed Pseudotime analysis with original seurat object as well to compare the results. 
    - Seurat Clusters: 10 (GE) and 14 (GTE)
    - Monocle Clusters: 11(GE) and 13 (GTE)
    - Created UMAPs with lower resolutions for ease of cluster cell classification: 0.2-0.6
- 2021-08-11:
    - Heatmaps of DEGs (resolution0.7)
    - RidgePlots of DEGs (resolution0.7)
- GBMSC_GSEA-aug31: GSEA results (.tsv summary file) was generated and analyzed in `neuroblastoma` notebook. 
    - GBM-subtypes Mapped to UMAP (resolution0.3)
    - Barplot of GSEA results to (resolution0.3) clusters


Dir: GBMSCTE 
Script: neuroblastoma_te.ipynb
Associated results in backup: 

- 2021-09-02: 
    - Seurat object with GBM-subtype labels, GTE object also contains Retrotransposon gene ratio
    - UMAP of Retrotransposon ratio + BarPlots of mean expression levels (%)
    - QQplots to test normality/ skewness. 
    - Bar plots of sample origin type of cells (GBM / NCS dataset) in each cluster and cluster 7 alone. 
    - Volcano plot of DEGs in cluster 7
- 2021-10-14:
    - Redid DEG analysis with Volcano plots and barplots of TE (organized by family,class and combined)and ZNFs (organized by group). 
    - Identified SVA, L1HS and ALu elements most significantly expressed in Cluster7
    - Analyzed DEGs between GBM and NCS cells within cluster 7. 
    - `neuroblastoma_te` output for cluster 7 retrotransposon analysis for human-specific Alu elements. 
- 2022-01-24:
    - Stored objects from Gene Ontology of cluster 7 GBM + GBMSC DEG from 2021-10-14. Made GO bar plot figures
   - 2022-01-25: Violin plot of teRatio (% transposon counts) by cluster with kruskal-wallis tes & wilconxon test to compare multiple groups and means between paired samples,respectively. 
   - 2022-01-30: Trajectory Analysis with GBM+GSC (hg38+retrotransposon mapped reads) with cluster6 (proneural cells) set as root cells. Making for BIN6005 poster presentation.
      - 2022-02-10: Using par1 object (cells ordered in pseudotime), the root cells were selected automatically given a the cluster (function finds the cells that come the earliest in the timeline). See `bin/r_gbmsc_monocle3_TEc6.R` for script.
- GBM_NSC_GSEA2022jan27: GSEA analysis on all clusters with senescence and cell type genesets.
- 2022-02-10: `bin/gbmsc_BIN6005posterfigs.R`
   - made figures without labels so that they can be adjusted in the powerpoint itself. UMAP seurat cluster (PC20r0.3), teRatio FeaturePlot, and violin plots. 
   - 2022-02-12: finished violin plot for retrotransposon expression; cleaned up Volcano Plot with gglot2; Started to clean up Gene Ontology (GO) terms in bar plot to pick only the most general terms related to cancer. 
- 2022-07-20: GE Clustered Gene analysis (csv, and plot). 

Dir: grant_compBrainDiseases
Script: grant_computationOfBrainDiseases_figs.ipynb 

- 2021-11-11:
    - Grant figures for Computational study of brain diseases. Contains variations of sample distribution, cluster PC20r0.3, GSEA subtypes, and Percentage of Retrotransposons expressed per cell. 
- 2022-05-18: 
    - `data-prep.ipynb` - preparing dataset for training the model.
    - `msm-model.ipynb` - creating and training the model.

Dir: GBM-GSC neuroblastoma TE
Script: FRQS_figs.R

- 2022-09-19: Dot plot drafts for FRQS poster. Describes GBM subtype marker differential expression across GSEA GBM subtype labeled clusters. 
- 2022-09-20: GO Term dot plot for cancer-related terms between TE-enriched cluster (7 in gte), GBM subset, and NSC subset (See `RedoGOanalysis_GBMvsNSC.ipynb`). Volcano plot of TE-enriched cluster (7). 

Dir: GBMGSCTE_dataprep
Scripts: data-prep.ipynb, data-prep.rmd, cnv_analysis.rmd

- 2022-06-29, 2022-07-04: Tried to merge the metadata of the original GBM dataset with my own to identify tumor cells, but they didn't match as indicated by Dr. Tanaka's table (see Slack). He ended up sending me the merged metadata instead of me figuring it out on my own. 


## Badhuri et al. 2020 Healthy Brain Controls (snRNAseq):

Dir: Brain_ctrl
Script: healthybrains.ipynb and healthybrain2.ipynb (UMAP testing)
Associated results in backup:

- ~~2021-06-15~~: Seurat object with quality control metrics
- 2021-06-23: ~~Seurat object with sample metadata~~ + Quality control figures
- 2021-06-30: Redid loading, quality control and filtering of cells. Filtered based on threshold over GTE dataset. 
- 2021-07-22: Seurat object with sample metadata and quality control metrics + QC figures from Seurrat vignette + QC figures after filtering, Top10 Variable Features after performing FindVariable(), UMAP figures at PC5-16 and resolution (0.2-1.2)
- 2021-08-05: 
    - Made new Seurat objects to match closer with the pre-processing and quality control performed in Bhaduri etal. 2020. 
    - Returning to Seurat object that was first made, saved ElbowPlots
    - UMAPs: chose PC14_r0.2 and r0.4
    - DEG: dataframe of markers for each cluster
- 2021-08-06: 
    - Barplots to determine cell types in each cluster 
    - Barplots showing sample composition per cluster
    - Barplots showing cluster composition for each sample 
    - Feature plots of brain tissue cell types. 
- 2021-08-11:
    - Heatmaps of brain markers & DEGs
    - RidgePlots of brain markers & DEGs

## Badhuri et al. 2020 Young Adult Healthy Brain Controls (snRNAseq):

Dir: BrainYA_ctrl
Script: yabrain.ipynb
Associated results in backup:

- 2021-08-31: ~~Seurat object~~(accidentally deleted, redone on 2021-11-25 using a seed) with quality control metrics, QC figures, FindVariable ScatterPlots, PCA QC, UMAP, Explore Known Brain Markers 
- 2021-09-03: DEG table of results, DEG Heatmaps, GO analysis (BarPlots)
- 2021-11-25: Redid progress on 2021-08-31, DEG and GO performed at a later date. 

## Merged Disease and Healthy Control Datasets

Dir: GBMGSC_BrainYA_merge
Script: gbmsc-yabrain_merged.ipynb
Associated results in backup:


