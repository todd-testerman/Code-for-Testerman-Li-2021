#!/bin/bash

filename='file_list.txt'

while read line; 
do
echo $line
metaphlan $line --bowtie2out $line.bowtie2.bz2 --nproc 16 --input_type fastq -o profiled_metagenome_$line.txt
done < $filename