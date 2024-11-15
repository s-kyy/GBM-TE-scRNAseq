#!/usr/bin/bash
vars=(./*/*_R*.fastq)
for DIR in "${vars[@]}"
do
    echo "zipping ${DIR}"
    gzip -v ${DIR} 
done
echo "FastQ compression complete"

vars=(./*/*_R*.fastq.gz)
for DIR in "${vars[@]}"
do
    echo "testing compression integrity ${DIR}"
    gzip -t ${DIR} 
done
echo "Compression integrity testing is complete"