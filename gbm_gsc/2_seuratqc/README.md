# Seurat Quality Control steps

Single cell analysis will be performed on two feature-cellbarcode matrices. One contains only the human reference genome annotations (**GE**), while the other also inlcudes retrotransposon mapped reads (**GETE**). 

### Setup

- Linux: CentOS7 ([StdEnv/2020](https://docs.alliancecan.ca/wiki/Standard_software_environments#StdEnv/2020))
- R 4.0.2
- renv 0.13.2
- BiocManager 1.30.12
- Seurat 4.0.1
- See comprehensive list of requirements in: `renv.lock`

Note: your setup should have enough memory to load the datasets into R and perform downstream proprocessing and analysis steps (16-32Gb )

### `1-createSeuratObj-GTE.R`

Script that aligns samples between human transcriptome-mapped reads to samples mapped to retrotransposon annotations.

One script was especially made for Bhaduri et al., 2020 GBM samples (`1-createSeuratObj-bhaduri.R`), because this data was analyzed before the automated script used for Wang et al., 2020 GBM samples and Bhaduri et al., 2020 Healthy samples was finalized.

In the **Load Datasets** section, replace `samples.csv`, `matrix.path` and `matrix.path.TE` with appropriate filepaths for the sampleIDs, as well as human transcriptome-mapped and human retrotransposon-mapped 10X count matrices. 

Values to evaluate the quality of the datasets is also computed in this script: 
- **Ratio of Genes (nGene) per UMI (nUMI) detected** - the greater the ratio, the more complex the dataset
- **Ratio of mitochondrial genes expressed** - the smaller the more likely a live cell at the time of dissociation (higher quality)

Commands and scripts used to create seurat objects for each dataset.

```bash
cd ./gbm_gsc/2_seuratqc

# 20210329_bhaduriGBM
Rscript --no-save --no-restore --verbose ./1-createSeuratObj-bhaduri.R >bhadurigbm_1.out 2>&1 

# 20210611_wangGBM
Rscript --no-save --no-restore --verbose ./1-createSeuratObj.R ../0_downloads/2021-05-28_wang/samples.csv \ 
../1_scrna-seq_mapping/2021-06-10_wang_aggr_ge \
../1_scrna-seq_mapping/2021-06-10_wang_aggr_te \
wangGBM >wanggbm_1.out 2>&1 

#20230320_healthy
Rscript --no-save --no-restore --verbose ./1-createSeuratObj.R ../0_downloads/2023-03-06_bhaduri_healthy/samples.csv \ 
../1_scrna-seq_mapping/2023-03-06_healthy_aggr_ge \
../1_scrna-seq_mapping/2023-03-06_healthy_aggr_te \
healthy >healthy_1.out 2>&1 
```
### `2-mergeSeuratObj.R`

Commands and scripts used to merge GBM datasets from Bhaduri et al., 2020 and Wang et al., 2020

```bash
cd ./gbm_gsc/2_seuratqc

# BhaduriGBM_GTE + WangGBM_GTE ./20230611_merged_bhaduri_GBM_SC/merged_GBM_GSC_gte.rds
Rscript --no-save --no-restore --verbose ./20210329_bhaduriGBM/seurat_obj/gte.rds ./20210611_wangGBM/seurat_obj/gte.rds GBM SC >merge_bhaduri_wang_gte.out 2>&1 

# BhaduriGBM_GE + WangGBM_GE => ./20230611_merge_bhaduri_GBM_SC/merged_GBM_GSC._ge.rds
Rscript --no-save --no-restore --verbose ./20210329_bhaduriGBM/seurat_obj/ge.rds ./20210611_wangGBM/seurat_obj/ge.rds GBM SC >merge_bhaduri_wang_ge.out 2>&1 
```

### `2-qcfigs.R`

Commands to generate quality control figures before quality control steps
```bash

```

Commands to generate quality control figures after filtering out low quality cells and doubets. 
```bash
```

### `3-normalize2umap.R`

```bash
```