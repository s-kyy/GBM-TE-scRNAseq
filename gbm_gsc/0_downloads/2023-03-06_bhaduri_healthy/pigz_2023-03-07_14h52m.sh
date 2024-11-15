#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=xxx
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --array=0-25

cd ./

if [[ $SLURM_ARRAY_TASK_ID == 0 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262920/SRR9262920_S9_L009_R1_001.fastq 
pigz -vt ./SRR9262920/SRR9262920_S9_L009_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 1 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262920/SRR9262920_S9_L009_R2_001.fastq 
pigz -vt ./SRR9262920/SRR9262920_S9_L009_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 2 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262922/SRR9262922_S10_L010_R1_001.fastq 
pigz -vt ./SRR9262922/SRR9262922_S10_L010_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 3 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262922/SRR9262922_S10_L010_R2_001.fastq 
pigz -vt ./SRR9262922/SRR9262922_S10_L010_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 4 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262923/SRR9262923_S11_L011_R1_001.fastq 
pigz -vt ./SRR9262923/SRR9262923_S11_L011_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 5 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262923/SRR9262923_S11_L011_R2_001.fastq 
pigz -vt ./SRR9262923/SRR9262923_S11_L011_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 6 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262932/SRR9262932_S12_L012_R1_001.fastq 
pigz -vt ./SRR9262932/SRR9262932_S12_L012_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 7 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262932/SRR9262932_S12_L012_R2_001.fastq 
pigz -vt ./SRR9262932/SRR9262932_S12_L012_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 8 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262937/SRR9262937_S13_L013_R1_001.fastq 
pigz -vt ./SRR9262937/SRR9262937_S13_L013_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 9 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262937/SRR9262937_S13_L013_R2_001.fastq 
pigz -vt ./SRR9262937/SRR9262937_S13_L013_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 10 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262938/SRR9262938_S5_L005_R1_001.fastq 
pigz -vt ./SRR9262938/SRR9262938_S5_L005_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 11 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262938/SRR9262938_S5_L005_R2_001.fastq 
pigz -vt ./SRR9262938/SRR9262938_S5_L005_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 12 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262946/SRR9262946_S8_L008_R1_001.fastq 
pigz -vt ./SRR9262946/SRR9262946_S8_L008_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 13 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262946/SRR9262946_S8_L008_R2_001.fastq 
pigz -vt ./SRR9262946/SRR9262946_S8_L008_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 14 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262956/SRR9262956_S6_L006_R1_001.fastq 
pigz -vt ./SRR9262956/SRR9262956_S6_L006_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 15 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9262956/SRR9262956_S6_L006_R2_001.fastq 
pigz -vt ./SRR9262956/SRR9262956_S6_L006_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 16 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9264382/SRR9264382_S1_L001_R1_001.fastq 
pigz -vt ./SRR9264382/SRR9264382_S1_L001_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 17 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9264382/SRR9264382_S1_L001_R2_001.fastq 
pigz -vt ./SRR9264382/SRR9264382_S1_L001_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 18 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9264383/SRR9264383_S3_L003_R1_001.fastq 
pigz -vt ./SRR9264383/SRR9264383_S3_L003_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 19 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9264383/SRR9264383_S3_L003_R2_001.fastq 
pigz -vt ./SRR9264383/SRR9264383_S3_L003_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 20 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9264385/SRR9264385_S7_L007_R1_001.fastq 
pigz -vt ./SRR9264385/SRR9264385_S7_L007_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 21 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9264385/SRR9264385_S7_L007_R2_001.fastq 
pigz -vt ./SRR9264385/SRR9264385_S7_L007_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 22 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9264388/SRR9264388_S4_L004_R1_001.fastq 
pigz -vt ./SRR9264388/SRR9264388_S4_L004_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 23 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9264388/SRR9264388_S4_L004_R2_001.fastq 
pigz -vt ./SRR9264388/SRR9264388_S4_L004_R2_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 24 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9264389/SRR9264389_S2_L002_R1_001.fastq 
pigz -vt ./SRR9264389/SRR9264389_S2_L002_R1_001.fastq.gz
fi

if [[ $SLURM_ARRAY_TASK_ID == 25 ]] ; then
pigz -8 -b 1024 -vk -p 8 ./SRR9264389/SRR9264389_S2_L002_R2_001.fastq 
pigz -vt ./SRR9264389/SRR9264389_S2_L002_R2_001.fastq.gz
fi

