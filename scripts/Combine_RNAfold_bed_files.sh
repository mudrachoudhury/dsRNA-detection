#!/bin/bash

#This script combines separate chromosomes together for each merge distance


declare -a distances=("1000")
#declare -a distances=("100" "150" "200" "400" "600" "1000")
for dist in "${distances[@]}"
do
	echo $dist
	input_files=$(ls '/home/mudrachoudhury/Xiaolab/dsRNA_detection_2020/data/Merged_'$dist'_EERs/Merged_'$dist'_split_chromosomes.chr'*'.output_w_RNA_fold.bed')
	output_file='/home/mudrachoudhury/Xiaolab/dsRNA_detection_2020/data/Merged_'$dist'_EERs/Merged_'$dist'_combined_output_w_RNA_fold.bed'
	rm $output_file
	cat $input_files >> $output_file
done
