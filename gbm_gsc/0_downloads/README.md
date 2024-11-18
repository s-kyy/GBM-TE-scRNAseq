# Download Datasets

## Download Bhaduri et al., 2020 GBM scRNA-seq 
Date: 2021-03-02

Install bamtofastq-1.3.2

```bash
cd ~/bin
wget https://cf.10xgenomics.com/misc/bamtofastq-1.3.2
chmod 700 bamtofastq-1.3.2

# ensure bin directory or executable directory is included in PATH (~/.bash_profile or ~/.bashrc)
export PATH=$PATH:~/bin

source ~/.bash_profile

# test
bamtofastq-1.3.2 -h
```

Download files

```bash
# download
wget -i /bhaduri2020_links.txt

# run convertbamtofastq scripts (slurm job)
sbatch convertBAM2FASTQ_0-7_9.sh  
sbatch convertBAM2FASTQ_8_10-11.sh
```

NOTE: The output directories (ie. data_##) should not exist prior to running convertBAM2FASTQ scripts.

## Download Wang et al., 2020 GBM CD133+/+ scRNA-seq
Date: 2021-05-28

Tools: sra-toolkit/2.10.8

```bash
# load tool on Cedar (Digital Research Alliance of Canada)
module load StdEnv/2020
module load gcc/9.3.0
module load sra-toolkit/2.10.8

chmod +x wang2020_download.sh
./wang2020_download.sh &> wang2020.out
```

## Download Bhaduri et al., 2020 Healthy snRNA-seq
Date: 2023-03-06

Tools: sra-toolkit/2.10.8, pigz v2.7

```bash
# load tool on Cedar (Digital Research Alliance of Canada)
module load StdEnv/2020
module load gcc/9.3.0
module load sra-toolkit/2.10.8

chmod +x bhaduri2020_healthy_download.sh
chmod +x bhaduri2020_healthy_zip.sh
cd ./2023-03-06_bhaduri_healthy
./bhaduri2020_healthy
./bhaduri2020_healthy_download_rename.sh ./samples.csv > tmp.out
./make_pigz_job.py
sbatch ./pigz_2023-03-07_14h52m.sh

# usage: make_pigz_job.py [-h] -i INPUT_FOLDER [-c CORES] [-m MEM]

# Produce script for cellranger count function

# options:
#   -h, --help       show this help message and exit
#   -i INPUT_FOLDER  fastq filepath (default: None)
#   -c CORES         local cores value used in SBATCH (e.g. 8, 12, 16). larger the value the higher greater the faster
#                    the process (default: 8)
#   -m MEM           mempercore value used in SBATCH (e.g. 10, 12, 15) (default: 10)
```

Note: File compression takes around 30min to 4h depending on file size. The generated script used in this project is provided in this repo: `pigz_2023-03-07_14h52m.sh`. Absolute paths were modified to relative paths for privacy purposes. 

# Download reference annotations

## Download Human Transcriptome Reference Annotations (GRCh38)

From the cellranger v6 downloads pages, download the human reference genome GRch38 2020. 

```bash
mkdir ./gbm_gsc/0_downloads/hg38-ref
cd ./gbm_gsc/0_downloads/hg38-ref
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
```

or execute script: `refdata_download.sh` (also includes transposable element annotations)

Install the cellranger v3.0.2 software: [https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0/). 
Ensure path to installed software is exported to PATH environment variable.

## Download and Make Retrotransposon Reference Annotations

Downloaded GTF annotations for transposable elements extracted with RepeatMasker from the human Ensembl genome assembly GRCh38 available from [Dr. Molly Gale Hammell's Lab](https://hammelllab.labsites.cshl.edu/software/#TEtranscripts).

```bash
mkdir ./gbm_gsc/0_downloads/te-ref
cd ./gbm_gsc/0_downloads/te-ref
wget https://www.dropbox.com/sh/1ppg2e0fbc64bqw/AAC5dKo3gsv9j-6eq5DHFJdha/GRCh38_Ensembl_rmsk_TE.gtf.gz?dl=1
mv GRCh38_Ensembl_rmsk_TE.gtf.gz?dl=1 GRCh38_Ensembl_rmsk_TE.gtf.gz
gzip -d GRCh38_Ensembl_rmsk_TE.gtf.gz 
```

or execute script: `refdata_download.sh` (also includes human reference annotations)

Next, in the `/0_downloads` run the following two scripts:

- `te-ref_filter.sh`: Filtered transposable element transcripts uniquely in the GTF file and not in the human reference `fa.fai` file. 
- `te-ref_mkref_job.sh`: Run mkref from cellranger v3.0.2
