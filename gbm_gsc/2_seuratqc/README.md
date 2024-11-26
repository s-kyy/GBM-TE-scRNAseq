# Seurat Quality Control steps

Single cell analysis will be performed on two feature-cellbarcode matrices. One contains only the human reference genome annotations (**GE**), while the other also inlcudes retrotransposon mapped reads (**GETE**). 

### Setup

- Linux: CentOS7 ([StdEnv/2020](https://docs.alliancecan.ca/wiki/Standard_software_environments#StdEnv/2020))
- R 4.0.2
- renv 0.13.2
- BiocManager 1.30.12
- Seurat 4.0.1
- ggplot2_3.3.5
- scales_1.1.1
- See comprehensive list of requirements in: `renv.lock`

Note: your setup should have enough memory to load the datasets into R and perform downstream proprocessing and analysis steps (16-32Gb )

### `1-createSeuratObj-GTE.R`

Script that aligns samples between human transcriptome-mapped reads to samples mapped to retrotransposon annotations.

One script was especially made for Bhaduri et al., 2020 GBM samples (`1-createSeuratObj-bhaduri.R`), because this data was analyzed before the automated script used for Wang et al., 2020 GBM samples and Bhaduri et al., 2020 Healthy samples was finalized.

In the **Load Datasets** section, replace `samples.csv`, `matrix.path` and `matrix.path.TE` with appropriate filepaths for the sampleIDs, as well as human transcriptome-mapped and human retrotransposon-mapped 10X count matrices. 

Values to evaluate the quality of the datasets is also computed in this script: 
- **Ratio of Genes (nGene) per UMI (nUMI) detected** (gene novelty score) - the greater the ratio, the more complex the dataset
- **Ratio of mitochondrial genes expressed** - the smaller the more likely a live cell at the time of dissociation (higher quality)
- **Median Absolute Deviations (MADs)** -- To detect if barcodes contain more than one cell, we calculate the MADs for the number of genes and UMIs detected per cell. Those expressing >3-5 MADs in number of genes or UMIs are likely doublets ([You et al., 2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02552-3) and [Ocasio et al., 2019](https://www.nature.com/articles/s41467-019-13657-6)). An elevated MAD threhold is preferred for cancer datasets due to the presence of copy number variations and other biological sources that increase variability in gene expression levels. 

Commands and scripts used to create seurat objects for each dataset.

```bash
cd ./gbm_gsc/2_seuratqc

Rscript --vanilla ./1-createSeuratObj-bhaduri.R >bhadurigbm_1.out 2>&1 
  # ./gbm_gsc/2_seuratqc/20210329_bhaduriGBM/gte.rds
  # ./gbm_gsc/2_seuratqc/20210329_bhaduriGBM/ge.rds

Rscript --vanilla ./1-createSeuratObj.R ../0_downloads/2021-05-28_wang/samples.csv \ 
../1_scrna-seq_mapping/2021-06-10_wang_aggr_ge \
../1_scrna-seq_mapping/2021-06-10_wang_aggr_te \
wangGBM >wanggbm_1.out 2>&1 
  # ./gbm_gsc/2_seuratqc/20210611_wangGBM/wangGBM_gte.rds
  # ./gbm_gsc/2_seuratqc/20210611_wangGBM/wangGBM_ge.rds

Rscript --vanilla ./1-createSeuratObj.R ../0_downloads/2023-03-06_bhaduri_healthy/samples.csv \ 
../1_scrna-seq_mapping/2023-03-06_healthy_aggr_ge \
../1_scrna-seq_mapping/2023-03-06_healthy_aggr_te \
healthy >healthy_1.out 2>&1 
  # ./gbm_gsc/2_seuratqc/20230320_healthy/healthy_gte.rds
  # ./gbm_gsc/2_seuratqc/20230320_healthy/healthy_ge.rds
```

Example of Running script on Windows (Powershell 7)

```bash
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla 1-createSeuratObj.R "..\0_downloads\2021-05-28_wang\samples.csv" "..\1_scrna-seq_mapping\2021-06-10_wang_aggr_ge" "..\1_scrna-seq_mapping\2021-06-10_wang_aggr_te" "wangGBM" *>"wanggbm_1.out" 
```

### `2-mergeSeuratObj.R`

Commands and scripts used to merge GBM datasets from Bhaduri et al., 2020 and Wang et al., 2020

```bash
cd ./gbm_gsc/2_seuratqc

Rscript --vanilla 2-mergeSeuratObj.R \
./20210329_bhaduriGBM/gte.rds \ 
./20210611_wangGBM/gte.rds bhaduriGBM wangGBM >merge_bhaduri_wang_gte.out 2>&1 
  # ./20230611_merged_bhaduriGBM_wangGBM/merged_bhaduriGBM_wangGBM_gte.rds

Rscript --vanilla 2-mergeSeuratObj.R \
./20210329_bhaduriGBM/ge.rds \
./20210611_wangGBM/ge.rds bhaduriGBM wangGBM >merge_bhaduri_wang_ge.out 2>&1 
  # ./20230611_merge_bhaduriGBM_wangGBM/merged_bhaduriGBM_wangGBM_ge.rds
```

### `qcfigs.R`

This script generates the following figures: 
- Bar plot of cell count per sample
- Density plot of genes detected per cell by sample
- Violin plot of genes detected per cell by sample 
- Scatter plot of genes correlated with UMIs per cell by sample
- Density plot of expressed mitochondrial genes per cell by sample
- Density plot of genes per UMI ratio per cell by sample

Commands to generate quality control figures before quality control steps. 

```bash
cd ./gbm_gsc/2_seuratqc

Rscript --vanilla qc-figs.R \
./20230611_merged_bhaduriGBM_wangGBM/merged_bhaduriGBM_wangGBM_gte.rds >qcfigs_merged.out 2>&1 

Rscript --vanilla qc-figs.R \
./20230611_merge_bhaduriGBM_wangGBM/merged_bhaduriGBM_wangGBM_ge.rds >qcfigs_merged.out 2>&1 

Rscript --vanilla qc-figs.R \
./20230320_healthy/healthy_gte.rds >qcfigs_healthygte.out 2>&1 

Rscript --vanilla qc-figs.R \
./20230320_healthy/healthy_ge.rds >qcfigs_healthyge.out 2>&1 
```

Minor adjustments to output figure sizes were done between runs. 

### `3-normalize2umap.R`

Accepting two filepaths to a human transcriptome-mapped seurat object and combined transcriptome and retrotransposon-mapped seurat object as inputs, this script outputs four RDS files in the same folder after filtering low quality cells based on quality control measures from previous script (`xxx_qc.rds`), after integrating samples via FindIntegrationAnchors and IntegrateData Seurat functions(`xxx_integrated.rds`), after dimensionality reduction with PCA (`xxx_integrated.umap.rds`) and after finding clusters through FindNeighbors and FindClusters Seurat functions at resolution 0.3 and 0.4 (`xxx_integrated.umap.clsutered.rds`). 

The following figures are also generated:
- `xxx_variablegenes.tiff` - Scatter plot individual genes described by average expression and standardized variance. 
- `xxx_UMAP-sample.tiff` - UMAP vizualization with cells annotated by sample
- `xxx_UMAP-cluster.tiff` - UMAP vizualization with cells annotated by cluster

```bash
cd ./gbm_gsc/2_seuratqc

Rscript --vanilla 2-mergeSeuratObj.R \
./20230611_merged_bhaduriGBM_wangGBM/merged_bhaduriGBM_wangGBM_gte.rds \ 
./20230611_merge_bhaduriGBM_wangGBM/merged_bhaduriGBM_wangGBM_ge.rds >normalize_bhaduri_wang.out 2>&1 
# 20230611_merged_bhaduriGBM_wangGBM/ -> location of saved RDS files
# 20230611_merged_bhaduriGBM_wangGBM/figs/ -> location of figures 

Rscript --vanilla 2-mergeSeuratObj.R \
./20230320_healthy/healthy_gte.rds \ 
./20230320_healthy/healthy_ge.rds >normalize_healthy.out 2>&1 
# 20230320_healthy/ -> location of saved RDS files
# 20230320_healthy/figs/ -> location of figures
```

<!-- Commands to generate quality control figures after filtering out low quality cells and doubets. 
```bash

``` -->
