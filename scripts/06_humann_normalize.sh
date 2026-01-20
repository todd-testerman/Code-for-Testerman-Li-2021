#!/bin/bash

filename='file_list_humann_normalize.txt'

while read line; 
do
echo $line
humann_renorm_table --input "$line"_\genefamilies.tsv --output "$line"_\genefamilies_relab.tsv --units relab --update-snames
done < $filename