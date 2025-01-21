# Seurat Clustering and Exploratory Analysis

### Requirements

- viridis 0.6.5
- future 1.34.0
- future.apply 1.11.3
- inferCNV 
- monocle3 1.3.1 with devtools. (see section below)
- SeuratWrappers 0.3.0. requires...
  - ellipsis_0.3.2
  - pillar_1.6.2
  - viridisLite 0.4.0
  - vctrs_0.3.8
<!-- - ComplexHeatmap 2.14 -->

#### Installing monocle3 on RHEL/CENTOS7 cluster environment

Load software

```bash
# on Niagara system
module load CCEnv arch/avx2 StdEnv/2020 gcc/9.3.0 gdal/3.5.1 geos/3.10.2 r/4.0.2
# on Cedar system
module load StdEnv/2020 gcc/9.3.0 gdal/3.5.1 geos/3.10.2 r/4.0.2
```

Install packages in R

```R
install.packages("devtools")
install.packages("remotes")
remotes::install_version("BiocManager", version="1.30.12")
library(BiocManager)
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor'))
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)

devtools::install_github('cole-trapnell-lab/leidenbase') #v0.1.9
remotes::install_version("matrixStats", version="0.58.0")

### Install sf
### ensure gdalv3.5.1+, geos/3.10.2, udunits/2.2.28 are available (and ideally no other version available)
###  Installation when multiple gdal versions are installed: https://github.com/r-spatial/sf/issues/844#issuecomment-653935662
remotes::install_version("udunits2", version="0.13.2.1")
remotes::install_version("units", version="0.8-5")
remotes::install_version("sf", version = '1.0-19', args="--no-test-load") 
  # last arguement is a quick fix for a bug in R that fails installation of sf library when multiple gdal versions are installed
  # Prior version: remotes::install_version("sf", version="0.3.2")

### Install spdep
remotes::install_version("spdep", version="1.1-8") 

### Install monocle3
devtools::install_github("cole-trapnell-lab/monocle3", ref="1.0.0")
  # OR 
  # devtools::install_github('cole-trapnell-lab/monocle3') 
  # devtools::install_github("cole-trapnell-lab/monocle3", ref="v1.3.1")

### Close and re-open new R session

### Install SeuratWrappers
remotes::install_version("R.utils", version="2.11.0")
devtools::install_github('satijalab/seurat-wrappers', ref="8510069") # v0.3.0 # still need to do (2025-01-15)
  # https://github.com/satijalab/seurat-wrappers/commit/8510069e76ae8f39c91250d67784e4b8ab4a9386
  # OR
  # devtools::install_github('satijalab/seurat-wrappers')
  # remotes::install_github('satijalab/seurat-wrappers@community-vignette') # v0.2.0
``` 

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

- UMAPs grouped or split by sample at two specified resolutions (default 0.3, 0.4)
- UMAPs grouped by cluster at two specified resolutions (default 0.3, 0.4)
- UMAP highlighting cells from female samples
- Ridge plot with cell cycle markers (PCNA, TOP2A, MCM6, MK67)
- [MA plots](https://www.jmp.com/support/downloads/JMPG101_documentation/Content/JMPGUserGuide/GR_G_0020.htm) of average expression of significantly expressed genes per cluster across log2 fold change and coloured by Bonferroni adjusted p values. 
- 

Commands to generate quality control figures after filtering out low quality cells. Also run again after each time object was integrated, and pre-processed to the clustering step (.e.g. removal of doublets, dead cells in gbm objects, and subsetting specific clusters, etc).

```bash
# Produce figures and summary CSVs after filtering low-quality cells. 
Rscript --vanilla ./2-clusteranalysis.R ./20241203_gbm_merged_ge_qc_integrated_umap/merged_ge_qc_integrated_umap_clustered.rds ./20241203_gbm_merged_ge_qc_integrated_umap/merged_ge_qc_integrated_umap_markers_0.3.csv ./20241203_gbm_merged_ge_qc_integrated_umap/merged_ge_qc_integrated_umap_markers_0.4.csv >gbm_ge_umap_figs.out 2>&1

Rscript --vanilla ./2-clusteranalysis.R ./20241203_gbm_merged_gte_qc_integrated_umap/merged_gte_qc_integrated_umap_clustered.rds ./20241203_gbm_merged_gte_qc_integrated_umap/merged_gte_qc_integrated_umap_markers_0.3.csv ./20241203_gbm_merged_gte_qc_integrated_umap/merged_gte_qc_integrated_umap_markers_0.4.csv >gbm_gte_umap_figs.out 2>&1

Rscript --vanilla ./2-clusteranalysis.R ./20241203_healthy_ge_qc_integrated_umap/ge_qc_integrated_umap_clustered.rds ./20241203_healthy_ge_qc_integrated_umap/ge_qc_integrated_umap_markers_0.3.csv ./20241203_healthy_ge_qc_integrated_umap/ge_qc_integrated_umap_markers_0.4.csv >healthy_ge_umap_figs.out 2>&1

Rscript --vanilla ./2-clusteranalysis.R ./20241203_healthy_gte_qc_integrated_umap/gte_qc_integrated_umap_clustered.rds ./20241203_healthy_gte_qc_integrated_umap/gte_qc_integrated_umap_markers_0.3.csv ./20241203_healthy_gte_qc_integrated_umap/gte_qc_integrated_umap_markers_0.4.csv >healthy_gte_umap_figs.out 2>&1
```

### `3-celltypeanalysis.R` (Preliminary)

Generate UMAPs colouring cells based on expression level of known celltype markers per cluster to identify cell type annotations. 

```bash
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3-celltypeanalysis.R ".\20241203_gbm_merged_ge_qc_integrated_umap\merged_ge_qc_integrated_umap_clustered.rds" ".\20241203_gbm_merged_ge_qc_integrated_umap\merged_ge_qc_integrated_umap_markers_0.3.csv" ".\20241203_gbm_merged_ge_qc_integrated_umap\merged_ge_qc_integrated_umap_markers_0.4.csv" *>gbm_ge_umap_celltypes_figs_dec03.out
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3-celltypeanalysis.R ".\20241203_gbm_merged_gte_qc_integrated_umap\merged_gte_qc_integrated_umap_clustered.rds" ".\20241203_gbm_merged_gte_qc_integrated_umap\merged_gte_qc_integrated_umap_markers_0.3.csv" ".\20241203_gbm_merged_gte_qc_integrated_umap\merged_gte_qc_integrated_umap_markers_0.4.csv" *>gbm_gte_umap_celltypes_figs_dec03.out
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3-celltypeanalysis.R ".\20241203_healthy_ge_qc_integrated_umap\ge_qc_integrated_umap_clustered.rds" ".\20241203_healthy_ge_qc_integrated_umap\ge_qc_integrated_umap_markers_0.3.csv" ".\20241203_healthy_ge_qc_integrated_umap\ge_qc_integrated_umap_markers_0.4.csv" *>healthy_ge_umap_celltypes_figs_dec03.out
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3-celltypeanalysis.R ".\20241203_healthy_gte_qc_integrated_umap\gte_qc_integrated_umap_clustered.rds" ".\20241203_healthy_gte_qc_integrated_umap\gte_qc_integrated_umap_markers_0.3.csv" ".\20241203_healthy_gte_qc_integrated_umap\gte_qc_integrated_umap_markers_0.4.csv" *>healthy_gte_umap_celltypes_figs_dec03.out
```

### Filter RBCs and Doublets `3.2-finddoublets.R` & `3.3-filterdoublets.R`

Red Blood cell gene expression was identified in step 3. 

Since RBCs are considered a contaminant in brain samples and healthy brain samples, we filtered these cells out by subsetting cells that expressed any classical RBC related genes [(Zhong et al., 2020)](https://www.nature.com/articles/s41586-019-1917-5).

Script removes RBCs, then identifies and removes possible doublets in each seurat object. We used DoubletFinder by McGinnis et al., (2019) due to its improved accuracy -at the cost of performance- compared to other available doublet annotation protocols in a recent comparitive review (Xi & Li, 2021).

Multiplet rates were estimated from (10X Genomics, 2019). We expect the following multiplet rates for each sample. 

- **GBM_Bhaduri2020** dataset loading rate was not mentioned, but their target capture rate was 2000 cells/sample, resulting in **0.8% multiplet rate.** 
- **GBM_Wang2020** dataset was loaded to 500 cells/uL at 21uL per sample, resulting in 10,500 cells/sample and **~2.6% multiplet rate**. 
- **Healthy_Bhaduri2020** dataset was loaded at 2000 nuclei/uL with a target capture rate of 3000 cells/sample, resulting in **~1.0% multiplet rate**.

Reference: 
- 10X Genomics. (2019). Single Cell 3’ Reagent Kits v2 User Guide RevF (Technical CG00052; p. 6). 10X Genomics. https://assets.ctfassets.net/an68im79xiti/RT8DYoZzhDJRBMrJCmVxl/6a0ed8015d89bf9602128a4c9f8962c8/CG00052_SingleCell3_ReagentKitv2UserGuide_RevF.pdf

```R
library(remotes) 
remotes::install_github(repo='chris-mcginnis-ucsf/DoubletFinder', ref = remotes::github_pull(176))
```

Once doublets are annotated and completed, `3.3-filterdoublets.R` script can be run to remove annotated doublets and repreprocess. 

#### (Test) Filtering of Dead Cells `3.5-filterDeadCells_gbm.R` with `3.5-filterDeadCells_gbm.sh`

Cluster 12 of GBM_GE object had elevated expression level of MALAT1, mitochondrial genes. We suspected it was a cluster of dead or hypoxic cells. They were filtered out with `3.5-filterDeadCells_gbm.R` with `3.5-filterDeadCells_gbm.sh` (not provided)  then re-integrated and reprocessed with `3.6-reintegrate-recluster_nia.R` and `3.6-re-int_niagara_*.sh` with the Seurat workflow as in [RBC and doublet filtering steps](#filter-rbcs-and-doublets-32-finddoubletsr--33-filterdoubletsr). 

- Genes expressed in less than 10 individual cells then re-integrated and prepocessed for further analysis. 
- Scaled data matrix was exported in a seperate object. 
- A higher resolution was used in FindClusters method (i.e. 0.5 and 0.6). 
- Alternatively, to run on cedar use scripts: `3.6-reintegrate_recluster.R` with `3.6-re-int_cedar.sh`

Results: clustering was less defined and a new hypoxic cluster was generated, suggesting that this cluster is inherent to the tumour dataset. Filter of dead cells was not retained for further analysis. 

#### Removal of lowly expressed genes `3.6-reintegrate-recluster_nia.R` and `3.6-re-int_niagara_*.sh`

- Genes expressed in less than 10 individual cells then re-integrated and prepocessed for further analysis. 
- Scaled data matrix was exported in a seperate object. 
- A higher resolution was used in FindClusters method (i.e. 0.5 and 0.6). 
- Alternatively, to run on cedar use scripts: `3.6-reintegrate_recluster.R` with `3.6-re-int_cedar.sh`

### Cell Type Annotation 

Run `2-clusteranalysis.R` and `3.4-finddoublets_validate.R` to generate summary tables (by sample or cluster) and figures. 

- UMAPs colouring cells by known cell markers
- DotPlots of cells against known marker gene expression and exports CSVs for number of cells per sample and per cluster to generate corresponding barplots. 

Then run `3.7-celltypeannotation.R` with `3.7-celltypeannotation.sh` to annotate clusters by cell type and create corresponding figures. 

```bash
## cluster analysis - post QC, filtRBC, filtDF (resolution 0.5, 0.6)
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\2-clusteranalysis.R ".\20250117_gbm_ge_filtDf_cluster\filtDf_cluster_int.rds" ".\20250117_gbm_ge_filtDf_cluster\filtDf_cluster_markers_0.5.csv" ".\20250117_gbm_ge_filtDf_cluster\filtDf_cluster_markers_0.6.csv" figs_clusteranalysis_ge_logFCall integrated_snn_res.0.5 integrated_snn_res.0.6 *>gbm_ge_filtDC_2figs0506_jan20.out

& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\2-clusteranalysis.R ".\20250117_gbm_gte_filtDf_cluster\filtDf_cluster_int.rds" ".\20250117_gbm_gte_filtDf_cluster\filtDf_cluster_markers_0.5.csv" ".\20250117_gbm_gte_filtDf_cluster\filtDf_cluster_markers_0.6.csv" figs_clusteranalysis_gte_logFCall integrated_snn_res.0.5 integrated_snn_res.0.6 *>gbm_gte_filtDC_2figs0506_jan20.out

& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\2-clusteranalysis.R ".\20250115_healthy_ge_filtDf_cluster\filtDf_cluster_int.rds" ".\20250115_healthy_ge_filtDf_cluster\filtDf_cluster_markers_0.5.csv" ".\20250115_healthy_ge_filtDf_cluster\filtDf_cluster_markers_0.6.csv" figs_clusteranalysis_ge integrated_snn_res.0.5 integrated_snn_res.0.6 *>healthy_ge_filtDC_2figs0506_jan15.out

& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\2-clusteranalysis.R ".\20250115_healthy_gte_filtDf_cluster\filtDf_cluster_int.rds" ".\20250115_healthy_gte_filtDf_cluster\filtDf_cluster_markers_0.5.csv" ".\20250115_healthy_gte_filtDf_cluster\filtDf_cluster_markers_0.6.csv" figs_clusteranalysis_gte integrated_snn_res.0.5 integrated_snn_res.0.6 *>healthy_gte_filtDC_2figs0506_jan15.out

## cell type validation - post QC, filtRBC, filtDF (resolution 0.5, 0.6)
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3.4-finddoublets_validate.R ".\20250117_gbm_ge_filtDf_cluster\filtDf_cluster_int.rds" figs_validatecelltypes_ge_res05 integrated_snn_res.0.5 *>gbm_ge_filtDC_34figs05_jan20.out
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3.4-finddoublets_validate.R ".\20250117_gbm_ge_filtDf_cluster\filtDf_cluster_int.rds" figs_validatecelltypes_ge_res06 integrated_snn_res.0.6 *>gbm_ge_filtDC_34figs06_jan20.out

& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3.4-finddoublets_validate.R ".\20250117_gbm_gte_filtDf_cluster\filtDf_cluster_int.rds" figs_validatecelltypes_gte_res05 integrated_snn_res.0.5 *>gbm_gte_filtDC_34figs05_jan20.out
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3.4-finddoublets_validate.R ".\20250117_gbm_gte_filtDf_cluster\filtDf_cluster_int.rds" figs_validatecelltypes_gte_res06 integrated_snn_res.0.6 *>gbm_gte_filtDC_34figs06_jan20.out

& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3.4-finddoublets_validate.R ".\20250115_healthy_ge_filtDf_cluster\filtDf_cluster_int.rds" figs_validatecelltypes_ge_res05 integrated_snn_res.0.5 *>healthy_ge_filtDC_34figs05_jan16.out
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3.4-finddoublets_validate.R ".\20250115_healthy_ge_filtDf_cluster\filtDf_cluster_int.rds" figs_validatecelltypes_ge_res06 integrated_snn_res.0.6 *>healthy_ge_filtDC_34figs06_jan16.out

& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3.4-finddoublets_validate.R ".\20250115_healthy_gte_filtDf_cluster\filtDf_cluster_int.rds" figs_validatecelltypes_gte_res05 integrated_snn_res.0.5 *>healthy_gte_filtDC_34figs05_jan16.out
& 'C:\Program Files\R\R-4.0.2\bin\Rscript.exe' --vanilla .\3.4-finddoublets_validate.R ".\20250115_healthy_gte_filtDf_cluster\filtDf_cluster_int.rds" figs_validatecelltypes_gte_res06 integrated_snn_res.0.6 *>healthy_gte_filtDC_34figs06_jan16.out
```

`3.6-celltype annotations` 


### `-cnvanalysis.R`

Install inferCNV with JAGS as a requirement. 

```R
library(BiocManager) # v3.14
BiocManager::install("infercnv") #v1.10.1
```

Environment:
```bash
module load StdEnv/2020
module load jags/4.3.2
module load r/4.1.0
```

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