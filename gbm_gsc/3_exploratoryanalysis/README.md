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

Generates figures for known celltype marker expression per cluster to identify cell type annotations. 

### `3.1-filterRBC_cluster.R`

Red Blood cell gene expression was identified in step 3. 

Since RBCs are considered a contaminant in brain samples and healthy brain samples, we filtered these cells out by subsetting cells that expressed any classical RBC related genes [(Zhong et al., 2020)](https://www.nature.com/articles/s41586-019-1917-5).

### `3.2-filterdoublets.R`

After removing uninformative cells and genes, we then removed possible doublets in each seurat object. 

We used DoubletFinder by McGinnis et al., (2019) due to its improved accuracy -at the cost of performance- compared to other available doublet annotation protocols in a recent comparitive review (Xi & Li, 2021).

```R
library(remotes) 
remotes::install_github(repo='chris-mcginnis-ucsf/DoubletFinder', ref = remotes::github_pull(176))
```

Once doublets are annotated and completed, rerun `1-pca_cluster.R` script. 

Multiplet rates were estimated from (10X Genomics, 2019). We expect the following multiplet rates for each sample. 

- **GBM_Bhaduri2020** dataset loading rate was not mentioned, but their target capture rate was 2000 cells/sample, resulting in **0.8% multiplet rate.** 
- **GBM_Wang2020** dataset was loaded to 500 cells/uL at 21uL per sample, resulting in 10,500 cells/sample and **~2.6% multiplet rate**. 
- **Healthy_Bhaduri2020** dataset was loaded at 2000 nuclei/uL with a target capture rate of 3000 cells/sample, resulting in **~1.0% multiplet rate**.

Reference: 
- 10X Genomics. (2019). Single Cell 3’ Reagent Kits v2 User Guide RevF (Technical CG00052; p. 6). 10X Genomics. https://assets.ctfassets.net/an68im79xiti/RT8DYoZzhDJRBMrJCmVxl/6a0ed8015d89bf9602128a4c9f8962c8/CG00052_SingleCell3_ReagentKitv2UserGuide_RevF.pdf

### `4-cnvanalysis.R`

Install CONICSmat (standalone software) and related requirements ([Müller et al., 2018](https://academic.oup.com/bioinformatics/article/34/18/3217/4979546?login=false)). Add paths to required software in `CONICS.cfg`. 

https://github.com/Neurosurgery-Brain-Tumor-Center-DiazLab/CONICS/wiki/Tutorial---CONICSmat;---Dataset:-SmartSeq2-scRNA-seq-of-Oligodendroglioma

```
module load StdEnv/2020
module load python/3.10.2
module load perl/5.30.2
module load samtools/1.17
module load bedtools/2.30.0
```
```
beanplotv1.3.1 (2022-04-09) - https://www.rdocumentation.org/packages/beanplot/versions/1.3.1
- Miktex (pdflatex in PATH)
- qpdf (qpdf present in PATH)
- Ghostscript in PATH
- pandoc
```

```
path_to_python="python"
path_to_samtools="samtools"
path_to_bedtools="bedtools"
path_to_rscript="Rscript"
path_to_bamreadcount="bam-readcount"
mappingCutoff=40
readCutoff=50000
fdrCutoff=0.05
genome="/path/to/genome.fa"
```