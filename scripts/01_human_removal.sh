#!/bin/bash

input="sample_list.txt";

while IFS= read -r line 
do
	echo "$line";
	kneaddata --input Raw_Fastq/"$line"_R1.fastq.gz --input Raw_Fastq/"$line"_R2.fastq.gz -db /Volumes/iMacPro_Pegasus/knead_database --output kneaddata_output -t 2 --run-fastqc-start --run-fastqc-end --trimmomatic /Users/users/miniconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1
	mv kneaddata_output/"$line"_R1_kneaddata_paired_1.fastq Cleaned_FastQ/"$line"_R1_kneaddata_paired_1.fastq
	mv kneaddata_output/"$line"_R1_kneaddata_paired_2.fastq Cleaned_FastQ/"$line"_R1_kneaddata_paired_2.fastq
done < $input