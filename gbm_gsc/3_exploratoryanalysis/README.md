# Seurat Clustering and Exploratory Analysis

### `1-pca_cluster.R` 

Takes one path as input with an optional output sample name. Outputs a seurat object after cell clustering and two CSVs with unique markers from each cell cluster at two resolutions (0.3 and 0.4). 

The number of principal components selected was based on the cumulative percentage of variation greater than 90% and variation associated with the PC is less than 5%

The following figures were also generated:
- Elbow Plot of standard deviation and cumulative variation of calculated PCs.
- CSV table of the variation captured in each PCs, along with their cumulative variation. 
- UMAP plot of based on sample 

Commands used: 

```bash
cd gbm_gsc/3_exploratoryanalysis

Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/ge_qc_integrated_regressCC.rds >healthy_ge_qc_rcc_pca.out 2>&1
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/gte_qc_integrated_regressCC.rds >healthy_gte_qc_rcc_pca.out 2>&1
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/ge_qc_integrated_regressCCdiff.rds >healthy_ge_qc_rccdiff_pca.out 2>&1
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/gte_qc_integrated_regressCCdiff.rds >healthy_gte_qc_rccdiff_pca.out 2>&1
# ./202421128_healthy/ <-- Output folder

Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20210902_merged_qc_bhaduriGBM_wangGBM/merged_ge_qc_integrated_regressCC.rds >gbm_ge_qc_rcc_pca.out 2>&1
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20210902_merged_qc_bhaduriGBM_wangGBM/merged_gte_qc_integrated_regressCC.rds >gbm_ge_qc_rcc_pca.out 2>&1
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20210902_merged_qc_bhaduriGBM_wangGBM/merged_ge_qc_integrated_regressCCdiff.rds >gbm_ge_qc_rccdiff_pca.out 2>&1
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20210902_merged_qc_bhaduriGBM_wangGBM/merged_gte_qc_integrated_regressCCdiff.rds >gbm_gte_qc_rccdiff_pca.out 2>&1
# ./20210902_merged_qc_bhaduriGBM_wangGBM
```

### `2-clusteranalysis.R`

Figures Generated: 

- UMAPs by Sample & Cluster at both 0.3 and 0.4 resolutions
- UMAP highlighting cells from female samples

Commands to generate quality control figures after filtering out low quality cells and doubets. 

```bash

```