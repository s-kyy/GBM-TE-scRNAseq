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