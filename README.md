# GBM IDH-wt Subtype Transitioning via Single Cell RNA and Retrotransposon Expression 

Glioblastoma is the most prevalent brain tumor in aging male population with a very poor survival rate. 
Although prognostic molecular features have been identified over the paste decades 
few were clinically relevant to improve treatment options. A hypothesis that has been gaining 
traction in GBM research is the presence of cancer stem cells responsible for initiating, driving and maintaining brain lesions. 

This study seeks to better characterize GBM stem cells in relation to retrotranposon expression which 
have been shown to be expressed in active stem cells during embryonic development. 

Analysis, experimental and figure-generating code for my master's memoir: [https://github.com/s-kyy/msc-memoir](https://github.com/s-kyy/msc-memoir).

## Datasets

Single cell RNA-seq analysis was performed on a publicly available NCBI datasets.

Training sets:

- Bhaduri et al. 2020. Cell: Stem Cell [SRP227039](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP227039)
- Wang et al. 2020. Stem Cell Reports [SRP227166](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP227166)

Healthy Brain controls: 
- Bhaduri et al. 2020. Cell: Stem Cell [SRP132816](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP132816&o=acc_s%3Aa)

Validation sets:

TCGA sets:  

## Requirements

Mapping of single cell or single nuclei RNA-seq datasets was performed on linux-based cluster resources (CentOS7). Analysis of read counts was performed on CentOS7 Linux-based remote cluster and Windows 10 operating systems. 

- bamtofastq v1.3.2 - https://cf.10xgenomics.com/misc/bamtofastq-1.3.2
- sra-toolkit v2.10.8 - https://github.com/ncbi/sra-tools/releases/tag/2.10.8
- cellranger v6.0.0 - https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0
- cellranger v3.0.2 (TE mapping) - https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/
- R v4.0.2 
    - Seurat v4.0
- R v4.4 
  - inferCNV
  - Seurat v4.4.0
  - SeuratObject v4.1.4
- Python v3.9


See the `README.md` under each numbered directory for more detailed notes on software requirements and script documentation. 

## Citations

See master's memoir: [https://github.com/s-kyy/msc-memoir](https://github.com/s-kyy/msc-memoir). 
