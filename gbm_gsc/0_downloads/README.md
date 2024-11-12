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
