# Download Bhaduri et al., 2020 scRNA-seq 
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

# Download Wang et al., 2020 scRNA-seq
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

# Download Human Transcriptome Reference Annotations (GRCh38)

From the cellranger v6 downloads pages, download the human reference genome GRch38 2020. 

```bash
mkdir ./gbm_gsc/0_downloads/hg38-ref
cd ./gbm_gsc/0_downloads/hg38-ref
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
```

# Download and Make Retrotransposon Reference Annotations

Downloaded GTF annotations for transposable elements extracted with RepeatMasker from the human Ensembl genome assembly GRCh38 available from [Dr. Molly Gale Hammell's Lab](https://hammelllab.labsites.cshl.edu/software/#TEtranscripts).

```bash
mkdir ./gbm_gsc/0_downloads/te-ref
cd ./gbm_gsc/0_downloads/te-ref
wget https://www.dropbox.com/sh/1ppg2e0fbc64bqw/AAC5dKo3gsv9j-6eq5DHFJdha/GRCh38_Ensembl_rmsk_TE.gtf.gz?dl=1
mv GRCh38_Ensembl_rmsk_TE.gtf.gz?dl=1 GRCh38_Ensembl_rmsk_TE.gtf.gz
gzip -d GRCh38_Ensembl_rmsk_TE.gtf.gz 
```

Filtered transposable element transcripts uniquely in the GTF file and not in the human reference `fa.fai` file. 

Used mkref function from cellranger v3.0.

## Remove the missing genes and save as a new GTF file, then run cellranger mkref
Filename: `slurm_cellrangermkref_RetroTEhg38.sh`

```shell
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-ytanaka
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G

cd /home/samkyy/projects/def-ytanaka/samkyy/resources/

grep -v -f uniq_missing GRCh38_Ensembl_rmsk_TE_v2.gtf > GRCh38_Ensembl_rmsk_TE_v3.gtf

/home/samkyy/bin/cellranger-6.0.0/cellranger mkref --genome=refdata_GRCh38-TE --fasta=/home/samkyy/projects/def-ytanaka/samkyy/resources/refdata-gex-GRCh38-2020-A/fasta/genome.fa
 --genes=/home/samkyy/projects/def-ytanaka/samkyy/resources/GRCh38_Ensembl_rmsk_TE_v23.gtf --nthreads=8 --memgb=10
```
- Takes a lot of memory: 10G --> 40G
