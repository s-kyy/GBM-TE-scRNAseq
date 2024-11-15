#!/usr/bin/bash

filename="$1"

echo ${pwd}

while read -r line || [[ -n "$line" ]] ; do
	runID="${line}"
	echo "prefetching ${runID}"

	prefetch ${runID}

	echo "fasterq-dump conversion to fastq files ${runID}"
	fasterq-dump --split-files ${runID} -O ./${runID} -t ./tmp &> ${runID}.out

	echo "fasterq-dump for ${runID} is complete"

done < "$filename"

#### RENAME ####
# > 30 year old samples
mv SRR9264382_1.fastq SRR9264382_S1_L001_R1_001.fastq
mv SRR9264382_2.fastq SRR9264382_S1_L001_R2_001.fastq

mv SRR9264389_1.fastq SRR9264389_S2_L002_R1_001.fastq
mv SRR9264389_2.fastq SRR9264389_S2_L002_R2_001.fastq

mv SRR9264388_1.fastq SRR9264388_S4_L004_R1_001.fastq
mv SRR9264388_2.fastq SRR9264388_S4_L004_R2_001.fastq

mv SRR9264383_1.fastq SRR9264383_S3_L003_R1_001.fastq
mv SRR9264383_2.fastq SRR9264383_S3_L003_R2_001.fastq

mv SRR9262938_1.fastq SRR9262938_S5_L005_R1_001.fastq
mv SRR9262938_2.fastq SRR9262938_S5_L005_R2_001.fastq

mv SRR9262956_1.fastq SRR9262956_S6_L006_R1_001.fastq
mv SRR9262956_2.fastq SRR9262956_S6_L006_R2_001.fastq

# <20 year old samples
mv SRR9264385_1.fastq SRR9264385_S7_L007_R1_001.fastq
mv SRR9264385_2.fastq SRR9264385_S7_L007_R2_001.fastq

mv SRR9262946_1.fastq SRR9262946_S8_L008_R1_001.fastq
mv SRR9262946_2.fastq SRR9262946_S8_L008_R2_001.fastq

mv SRR9262920_1.fastq SRR9262920_S9_L009_R1_001.fastq
mv SRR9262920_2.fastq SRR9262920_S9_L009_R2_001.fastq

mv SRR9262922_1.fastq SRR9262922_S10_L010_R1_001.fastq
mv SRR9262922_2.fastq SRR9262922_S10_L010_R2_001.fastq

mv SRR9262923_1.fastq SRR9262923_S11_L011_R1_001.fastq
mv SRR9262923_2.fastq SRR9262923_S11_L011_R2_001.fastq

mv SRR9262932_1.fastq SRR9262932_S12_L012_R1_001.fastq
mv SRR9262932_2.fastq SRR9262932_S12_L012_R2_001.fastq

mv SRR9262937_1.fastq SRR9262937_S13_L013_R1_001.fastq
mv SRR9262937_2.fastq SRR9262937_S13_L013_R2_001.fastq
