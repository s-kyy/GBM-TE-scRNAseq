# scRNA-seq Mapping FastQ reads to Reference Annotations

## Requirements

- cellranger v6.0.0 - [https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.0)
- cellranger v3.0.2 (TE mapping) - [https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/)
- Python 3.9.2 

Ensure path(s) to installed software is exported to PATH environment variable.

## Scripts

Overview of scripts in this directory (see sections below for details):

1. [Run `make_cellranger_job.py`](#1-run-make_cellranger_jobpy-)
2. [Run `./<date>_<dataset>_counts_<ref>/<dataset>_counts_<ref><data-time>.sh`](#2-run-date_dataset_counts_refdataset_counts_refdata-timesh-)
3. [Run `run_cellranger3.sh` + `run_cellranger6.sh`](#3-run-run_cellranger3sh--run_cellranger3sh-)
4. [Run cellranger aggregate scripts](#4-run-cellranger-aggregate-scripts-) 

### 1. Run `make_cellranger_job.py` [^](#scripts)

This python script generates a bash script which requests user-defined remote cluster resources (as SLURM commands), then calls  `run_cellranger#.sh` to run the count method using the defined version of cellranger. 

```bash
python ./make_cellranger_job.py 
# usage: make_cellranger_job.py [-h] -i INPUT_FOLDER -o OUT_FOLDER -r REF [-c CORES] [-m MEM] [-t TYPE]
# Produce script for cellranger count function

# options:
#   -h, --help       show this help message and exit
#   -i INPUT_FOLDER  fastq filepath (default: None)
#   -o OUT_FOLDER    output filepath for cellranger count (ie. "./") (default: None)
#   -r REF           reference annotations filepath (default: None)
#   -c CORES         local cores value used in SBATCH (e.g. 8, 12, 16). larger the value the higher greater the faster
#                    the process (default: 8)
#   -m MEM           mempercore value used in SBATCH (e.g. 10, 12, 15) (default: 10)
#   -t TYPE          cellranger version to use. (6 or 3) (default: 6)
#   -p PATTERN       SF... = 0 (default), SRR... = 1 (default: 0)
```

Note: 

List of commands

```bash
cd ./<date>_bhaduri_counts_ge
python ../make_cellranger_job.py -i ../../0_downloads/<date>_bhaduri -o ./ -r ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A -c 8 -m 10 -t 6 -p 0 &> tmp.out

cd ../<date>_bhaduri_counts_te
python ../make_cellranger_job.py -i ../../0_downloads/<date>_bhaduri -o ./ -r ../../0_downloads/te-ref/refdata_GRCh38-TE -c 8 -m 15 -t 3 -p 0 &> tmp.out

cd ../<date>_wang_counts_ge
python ../make_cellranger_job.py -i ../../0_downloads/<date>_wang -o ./ -r ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A -c 8 -m 15 -t 6 -p 1 &> tmp.out

cd ../<date>_wang_counts_te
python ../make_cellranger_job.py -i ../../0_downloads/<date>_wang -o ./ -r ../../0_downloads/te-ref/refdata_GRCh38-TE -c 8 -m 15 -t 3 -p 1 &> tmp.out

cd ../<date>_healthy_counts_ge
python ../make_cellranger_job.py -i ../../0_downloads/<date>_bhaduri_healthy -o ./ -r ../../0_downloads/hg38-refdata/refdata-gex-GRCh38-2020-A -c 8 -m 15 -t 6 -p 1 &> tmp.out

cd ../<date>_healthy_counts_te
python ../make_cellranger_job.py -i ../../0_downloads/<date>_bhaduri_healthy -o ./ -r ../../0_downloads/te-ref/refdata_GRCh38-TE -c 8 -m 15 -t 3 -p 1 &> tmp.out
```

*NOTE*: Generated scripts should use absolute paths. The scripts provided in this repository were modified to use relative paths for privacy. 

### 2. Run `./<date>_<dataset>_counts_<ref>/<dataset>_counts_<ref><data-time>.sh` [^](#scripts)

Output from python script #1. 

To run the generated output scripts, do so in within their output folders.

### 3. Run `run_cellranger3.sh` + `run_cellranger3.sh` [^](#scripts)

These bash scripts directly calls the appropriate version of cellranger using defined paramters from #2

### 4. Run cellranger aggregate scripts [^](#scripts)

Normalization was modified to none, and no secondary analysis was performed with cellranger since it will be done with Seurat. 

Run each script while in the `1_scrna-seq_mapping` directory. 

Resulting `<aggr_dir>/outs/count/feature_bc_matrix` will be used for subsequent preprocessing and analysis in R. 