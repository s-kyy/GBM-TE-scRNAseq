# scRNA-seq Mapping FastQ reads to Reference Annotations

## Requirements

- cellranger v6.0.0 - [https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0)
- cellranger v3.0.2 (TE mapping) - [https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/)
- Python 3.9.2 

Ensure path(s) to installed software is exported to PATH environment variable.

## Scripts

Overview of scripts in this directory (see sections below for details):

1. [`make_cellranger_job.py`](#1-run-make_cellranger_jobpy-)
2. [`./<date>_<dataset>_counts_<ref>/<dataset>_counts_<ref><data-time>.sh`](#2-run-date_dataset_counts_refdataset_counts_refdata-timesh-)
3. [`run_cellranger3.sh` + `run_cellranger6.sh`](#3-run_cellranger3sh--run_cellranger3sh)

### 1. Run `make_cellranger_job.py` [^](#scripts)

This python script generates a bash script which requests user-defined remote cluster resources (as SLURM commands), then calls  `run_cellranger#.sh` to run the count method using the defined version of cellranger. 

```bash
python ./make_cellranger_job.py 
# Usage: make_cellranger_job.py [fastqdir] [output_dir] [ref_genome] [localcores] [mempercore] [v3-v6] [output_name] 
```

Parameters
- `fastqdir` (path) = directory path containing list of fastq directories 
- `output_dir` (path) = directory path where cellranger count will output, path should already exist
- `ref_genome` (path) = directory path of the human reference annotations or transposable element reference annotations
- `localcores` (int) = total number of cores available (define less than maximum if performing locally)
- `mempercore` (int) = memory available per cpu (in Gb)
- `v3-v6` (int) = use cellranger 3 (map to transposon annotation) or 6 (human reference genome)
- `output_name` (str) = Prefix for job script file name. 

Example command:

```bash
cd ./2021-03-23_bhaduri_counts_ge
python ../make_cellranger_job.py ../../0_downloads/2021-03-02_bhaduri ./ ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A 8 15 6 bhaduri_counts_ge &> tmp.out
```

*NOTE*: Generated scripts will contain absolute paths. The generated scripts provided in this repository were modified to use relative paths for privacy. 



### 2. Run `./<date>_<dataset>_counts_<ref>/<dataset>_counts_<ref><data-time>.sh` [^](#scripts)

Output from python script #1. 

To run the generated output scripts, do so in within their output folders.

### 3. `run_cellranger3.sh` + `run_cellranger3.sh` [^](#scripts)

These bash scripts directly calls the appropriate version of cellranger using defined paramters from #2
