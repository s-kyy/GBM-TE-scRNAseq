# Seurat Quality Control steps

Single cell analysis will be performed on two feature-cellbarcode matrices. One contains only the human reference genome annotations (**GE**), while the other also inlcudes retrotransposon mapped reads (**GETE**). 

### Setup

- Linux: CentOS7 ([StdEnv/2020](https://docs.alliancecan.ca/wiki/Standard_software_environments#StdEnv/2020))
- R 4.0.2
- renv 0.13.2
- BiocManager 1.30.12
- Seurat 4.0.1
- See comprehensive list of requirements in: `renv.lock`

### `1-createSeuratObj-GTE.R`

Script that aligns samples between human transcriptome-mapped reads to samples mapped to retrotransposon annotations.

One script was especially made for Bhaduri et al., 2020 GBM samples (`1-createSeuratObj-bhaduri-GTE.R`), because this data was analyzed before the automated script used for Wang et al., 2020 GBM samples and Bhaduri et al., 2020 Healthy samples was finalized.

Values to evaluate the quality of the datasets is also computed in this script: 
- **Ratio of Genes (nGene) per UMI (nUMI) detected** - the greater the ratio, the more complex the dataset
- **Ratio of mitochondrial genes expressed** - the smaller the more likely a live cell at the time of dissociation (higher quality)

In the **Load Datasets** section, replace `samples.csv`, `matrix.path` and `matrix.path.TE` with appropriate filepaths for the sampleIDs, as well as human transcriptome-mapped and human retrotransposon-mapped 10X count matrices. 

### `2-qcfigs.R`



### `3-normalize2umap.R`