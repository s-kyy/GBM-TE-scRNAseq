# scRNA-seq Mapping FastQ reads to Reference Annotations

## Requirements

- cellranger v6.0.0 - [https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0)
- cellranger v3.0.2 (TE mapping) - [https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/)
- Python 3.9.2 

Ensure path(s) to installed software is exported to PATH environment variable.

## Scripts

Overview of scripts in this directory (see sections below for details):

1. [`make_cellranger_job.py`](#1-make_cellranger_jobpy)
2. [`run_cellranger3.sh` + `run_cellranger6.sh`](#2-run_cellranger3sh--run_cellranger3sh)

### 1. `make_cellranger_job.py` [^](#scripts)

This python script generates a bash script which requests user-defined remote cluster resources (as SLURM commands), then calls  `run_cellranger#.sh` to run the count method using the defined version of cellranger. 

```bash
python .\make_cellranger_job.py 
# Usage: make_cellranger_job.py [fastqdir] [output_dir] [ref_genome] [localcores] [mempercore] [v3-v6] [output_name] 
```

Example command:

```bash
cd ./2021-03-23_bhaduri_counts_ge
python ../make_cellranger_job.py ../../0_downloads/2021-03-02_bhaduri ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 6 bhaduri_counts_ge &> out
```

*NOTE*: Scripts used to run cellranger are provided in this repository. Absolute paths were modified for relative paths for privacy. 

### 2. `run_cellranger3.sh` + `run_cellranger3.sh` [^](#scripts)

This bash script directly calls the appropriate version of cellranger. 


