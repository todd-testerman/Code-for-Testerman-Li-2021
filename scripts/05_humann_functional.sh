#!/bin/bash

filename='file_list_humann.txt'

while read line; 
do
echo $line
humann --input concatenate_output/$line --output $line --threads 16
done < $filename