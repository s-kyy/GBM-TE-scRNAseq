# Seurat Clustering and Exploratory Analysis

### Requirements

- viridis 0.6.5
- viridisLite 0.4.1

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

Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/20241203_healthy_int_umap/ge_qc_integrated.rds >healthy_ge_qc_rcc_pca.out 2>&1
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20230320_healthy/20241203_healthy_int_umap/gte_qc_integrated.rds >healthy_gte_qc_rcc_pca.out 2>&1
# ./202421128_healthy/ <-- Output folder

Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20210902_merged_qc_bhaduriGBM_wangGBM/20241203_gbm_int_umap/merged_ge_qc_integrated.rds >gbm_ge_qc_rcc_pca.out 2>&1
Rscript --vanilla ./1-pca_cluster.R ../2_seuratqc/20210902_merged_qc_bhaduriGBM_wangGBM/20241203_gbm_int_umap/merged_gte_qc_integrated.rds >gbm_ge_qc_rcc_pca.out 2>&1
# ./20210902_merged_qc_bhaduriGBM_wangGBM
```

### `2-clusteranalysis.R`

Figures Generated: 

- UMAPs by Sample & Cluster at both 0.3 and 0.4 resolutions
- UMAP highlighting cells from female samples

Commands to generate quality control figures after filtering out low quality cells and doubets. 

```bash
Rscript --vanilla ./2-clusteranalysis.R ./20241203_gbm_merged_ge_qc_integrated_umap/merged_ge_qc_integrated_umap_clustered.rds ./20241203_gbm_merged_ge_qc_integrated_umap/merged_ge_qc_integrated_umap_markers_0.3.csv ./20241203_gbm_merged_ge_qc_integrated_umap/merged_ge_qc_integrated_umap_markers_0.4.csv >gbm_ge_umap_figs.out 2>&1

Rscript --vanilla ./2-clusteranalysis.R ./20241203_gbm_merged_gte_qc_integrated_umap/merged_gte_qc_integrated_umap_clustered.rds ./20241203_gbm_merged_gte_qc_integrated_umap/merged_gte_qc_integrated_umap_markers_0.3.csv ./20241203_gbm_merged_gte_qc_integrated_umap/merged_gte_qc_integrated_umap_markers_0.4.csv >gbm_gte_umap_figs.out 2>&1

Rscript --vanilla ./2-clusteranalysis.R ./20241203_healthy_ge_qc_integrated_umap/ge_qc_integrated_umap_clustered.rds ./20241203_healthy_ge_qc_integrated_umap/ge_qc_integrated_umap_markers_0.3.csv ./20241203_healthy_ge_qc_integrated_umap/ge_qc_integrated_umap_markers_0.4.csv >healthy_ge_umap_figs.out 2>&1

Rscript --vanilla ./2-clusteranalysis.R ./20241203_healthy_gte_qc_integrated_umap/gte_qc_integrated_umap_clustered.rds ./20241203_healthy_gte_qc_integrated_umap/gte_qc_integrated_umap_markers_0.3.csv ./20241203_healthy_gte_qc_integrated_umap/gte_qc_integrated_umap_markers_0.4.csv >healthy_gte_umap_figs.out 2>&1
```

### `3-celltypeanalysis.R`
### `3.1-filterRBC_cluster.R`
### `3.2-filterdoublets.R`

```R
library(remotes) 
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
```

Once doublets are annotated and completed, rerun `1-pca_cluster.R` script. 

### `4-cnvanalysis.R`

Install inferCNV v with JAGS as a requirement. 

```bash
module load StdEnv/2020
module load gcc/9.3.0
module load jags/4.3.2
module load r/4.1.0
```

Install inferCNV

```R
library(BiocManager) # v3.14
BiocManager::install("infercnv") #v1.10.1
```
