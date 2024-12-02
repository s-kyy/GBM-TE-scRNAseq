# Seurat Quality Control and Integration steps

Single cell analysis will be performed on two feature-cellbarcode matrices. One contains only the human reference genome annotations (**GE**), while the other also inlcudes retrotransposon mapped reads (**GTE**). 

### Setup

- Linux: CentOS7 ([StdEnv/2020](https://docs.alliancecan.ca/wiki/Standard_software_environments#StdEnv/2020))
- R 4.0.2
- renv 0.13.2
- BiocManager 1.30.12
- Seurat 4.0.1
- Matrix 1.3.-2
- tidyverse 1.3.1
- dplyr 1.07
- ggplot2 3.3.5
- ggpp 0.5.8-1
- ggpmisc 0.6.1
- scales 1.1.1
- conflicted 1.1.0
- fs 1.5.0
- See comprehensive list of requirements in: `renv.lock`

Note: your setup should have enough memory to load the datasets into R and perform downstream proprocessing and analysis steps (16-32Gb )

### `1-createSeuratObj.R`

Script that aligns samples between human transcriptome-mapped reads to samples mapped to retrotransposon annotations.

One script was especially made for Bhaduri et al., 2020 GBM samples (`1-createSeuratObj-bhaduri.R`), because this data was analyzed before the automated script used for Wang et al., 2020 GBM samples and Bhaduri et al., 2020 Healthy samples was finalized.

In the **Load Datasets** section, replace `samples.csv`, `matrix.path` and `matrix.path.TE` with appropriate filepaths for the sampleIDs, as well as human transcriptome-mapped and human retrotransposon-mapped 10X count matrices. 

Values to evaluate the quality of the datasets is also computed in this script: 
- **Ratio of Genes (nGene) per UMI (nUMI) detected** (gene novelty score) - the greater the ratio, the more complex the dataset
- **Ratio of mitochondrial genes expressed** - the smaller the more likely a live cell at the time of dissociation (higher quality)

Commands and scripts used to create seurat objects for each dataset.

```bash
cd ./gbm_gsc/2_seuratqc

Rscript --vanilla ./1-createSeuratObj-bhaduri.R >bhadurigbm_1.out 2>&1 
  # outputs:
  # ./gbm_gsc/2_seuratqc/20210329_bhaduriGBM/gte.rds
  # ./gbm_gsc/2_seuratqc/20210329_bhaduriGBM/ge.rds

Rscript --vanilla ./1-createSeuratObj.R ../0_downloads/2021-05-28_wang/samples.csv \ 
../1_scrna-seq_mapping/2021-06-10_wang_aggr_ge \
../1_scrna-seq_mapping/2021-06-10_wang_aggr_te \
wangGBM >wanggbm_1.out 2>&1 
  # outputs:
  # ./gbm_gsc/2_seuratqc/20210611_wangGBM/gte.rds
  # ./gbm_gsc/2_seuratqc/20210611_wangGBM/ge.rds

Rscript --vanilla ./1-createSeuratObj.R ../0_downloads/2023-03-06_bhaduri_healthy/samples.csv \ 
../1_scrna-seq_mapping/2023-03-06_healthy_aggr_ge \
../1_scrna-seq_mapping/2023-03-06_healthy_aggr_te \
healthy >healthy_1.out 2>&1 
  # outputs:
  # ./gbm_gsc/2_seuratqc/20230320_healthy/gte.rds
  # ./gbm_gsc/2_seuratqc/20230320_healthy/ge.rds
```

Example of Running script on Windows (Powershell 7)

```bash
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla 1-createSeuratObj.R "..\0_downloads\2021-05-28_wang\samples.csv" "..\1_scrna-seq_mapping\2021-06-10_wang_aggr_ge" "..\1_scrna-seq_mapping\2021-06-10_wang_aggr_te" "wangGBM" *>"wanggbm_1.out" 
```

### `2-mergeSeuratObj.R`

Commands and scripts used to merge GBM datasets from Bhaduri et al., 2020 and Wang et al., 2020. 

```bash
cd ./gbm_gsc/2_seuratqc

Rscript --vanilla 2-mergeSeuratObj.R \
./20210329_bhaduriGBM/gte.rds \ 
./20210611_wangGBM/gte.rds bhaduriGBM wangGBM >merge_bhaduri_wang_gte.out 2>&1 
  # outputs:
  # ./20230611_merged_bhaduriGBM_wangGBM/merged_gte.rds

Rscript --vanilla 2-mergeSeuratObj.R \
./20210329_bhaduriGBM/ge.rds \
./20210611_wangGBM/ge.rds bhaduriGBM wangGBM >merge_bhaduri_wang_ge.out 2>&1 
  # outputs:
  # ./20230611_merge_bhaduriGBM_wangGBM/merged_ge.rds
```

### `qcfigs_bhaduri.R` and `qcfigs_wang.R`

This script generates the following figures: 
- Bar plot of cell count per sample
- Density plot of UMI detected per cell by sample
- Density plot of genes detected per cell by sample
- Scatter plot of genes correlated with UMIs per cell by sample
- Density plot of expressed mitochondrial genes per cell by sample
- Density plot of genes per UMI ratio per cell by sample
- CSV of number of cells per sample with statistics on novelty score, as well as median absolute deviations of the number of genes detected, UMI detected and mitochondrial percentage. 

Commands to generate quality control figures before quality control steps. 

```bash
cd ./gbm_gsc/2_seuratqc

Rscript --vanilla qc-figs_bhaduriwang.R \
./20230611_merged_bhaduriGBM_wangGBM/merged_gte.rds >qcfigs_merged.out 2>&1 

Rscript --vanilla qc-figs_bhaduriwang.R \
./20230611_merge_bhaduriGBM_wangGBM/merged_ge.rds >qcfigs_merged.out 2>&1 

Rscript --vanilla qc-figs_healthy.R \
./20230320_healthy/gte.rds >qcfigs_healthygte.out 2>&1 

Rscript --vanilla qc-figs_healthy.R \
./20230320_healthy/ge.rds >qcfigs_healthyge.out 2>&1 
```

### `3-filter_gbm.R` and `3-filter_healthy.R`

Since the healthy samples are single nucleotide RNA-seq samples, slightly modified thresholds were used to filter low quality cells. We followed the thresholds used by Bhaduri et al., (2020) to filter snRNA-seq datasets. 

Some code has been commented out to filter possible doublets / triplets use 3-5 times the median absolute deviations of genes and UMI detected. However this was omitted in favor of using a unified R package to simulate and filter doublets (i.e. DoubletFinder). 

```bash
Rscript --vanilla 3-filter_gbm.R ./20230611_merge_bhaduriGBM_wangGBM/merged_ge.rds ./20230611_merge_bhaduriGBM_wangGBM/merged_gte.rds >filter_gbm.out 2>&1
  # outputs:
  # ./20230611_merge_bhaduriGBM_wangGBM/merged_ge_qc.rds
  # ./20230611_merge_bhaduriGBM_wangGBM/merged_gte_qc.rds

Rscript --vanilla 3-filter_gbm.R ./20230320_healthy/ge.rds ./20230320_healthy/gte.rds >filter_healthy.out 2>&1
  # outputs:
  # ./20230320_healthy/ge_qc.rds
  # ./20230320_healthy/gte_qc.rds
```

### `4-integrate.R`

This script accepts two filepaths to a seurat object (.rds) and list of cell cycle genes (.csv) as input while outputting two new seurat objects (.rds). One of which regresses out cell cycle related genes as to remove effects of mitosis and proliferation and ideally cluster based on more relevant cell types and tumours states. We also export intermediate seurat objects after computationally internsive steps (i.e. FindIntegrationAnchors and IntegrateData). 

Log normalization was used to standardize the raw counts. Top 2500 variable features were identified and used in FindIntegrationAnchors. When scaling and centering data S phase and G2/M phase scores were used to regress cell cycle effects. 

We included an example job script: `4-integrate_job.sh`. 

```bash
cd ./gbm_gsc/2_seuratqc

Rscript --vanilla ./4-integrate.R ./20230320_healthy/ge_qc.rds ./cellcycle_genes.csv >healthy_ge_qc_integrate.out 2>&1
Rscript --vanilla ./4-integrate.R ./20230320_healthy/gte_qc.rds ./cellcycle_genes.csv >healthy_gte_qc_integrate.out 2>&1
# 20230320_healthy/ --> location of saved RDS files
# 20230320_healthy/figs --> location of figures

Rscript --vanilla ./4-integrate.R ./20230611_merged_bhaduriGBM_wangGBM/merged_ge_qc.rds ./cellcycle_genes.csv >gbm_ge_qc_integrate.out 2>&1
Rscript --vanilla ./4-integrate.R ./20230611_merged_bhaduriGBM_wangGBM/merged_gte_qc.rds ./cellcycle_genes.csv >gbm_gte_qc_integrate.out 2>&1
# 20230611_merged_bhaduriGBM_wangGBM/ -> location of saved RDS files
# 20230611_merged_bhaduriGBM_wangGBM/figs/ -> location of figures 
```
