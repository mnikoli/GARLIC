#!/bin/bash

#
# Script to preprocess regulatory map like in the manuscript from GARLIC:
# Keeps only chromosomes 1-22,X,Y and removes lines with other chr ids, 
# Extends every region to a fixed size of 2000 nts.
# Place the script in the folder with your BED files and run it using ./regMapRegions_setup.sh.
# Make sure to change permission of the file so that it can be executed, using the command sudo chmod +x regMapRegions_setup.sh 
# 


for file in *.bed;
do

#echo $file
tmp_name=${file#*-}
tmp_name1=${tmp_name#*_}

new_name=${tmp_name1#*_} 
echo ${new_name%.bed}_input.csv

awk '$1!~/M/ && $1!~/GL/ && $1!~/_/ {print $1"\t"int(($3+$2)/2)"\t"int(($3+$2)/2)}' $file  | awk '{if ($2<=0) print $1"\t1\t"$3-$2; else print $0}' | awk ' {print $1"\t"$2-1000"\t"$3+1000}' | sortBed -i | mergeBed -i > "${new_name%.bed}_input.csv"

done

