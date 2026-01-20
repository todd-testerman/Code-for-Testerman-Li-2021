#!/bin/bash

filename='file_list_humann.txt'

while read line; 
do
echo $line
cat $line > "concatenate_output/merged_$line"
done < $filename