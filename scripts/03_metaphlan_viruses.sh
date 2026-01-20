#!/bin/bash

filename='file_list.txt'

while read line; 
do
echo $line
metaphlan $line.bowtie2.bz2 --nproc 16 --input_type bowtie2out --add_viruses -t rel_ab_w_read_stats -o abs_abun/profiled_metagenome_$line.txt
done < $filename