#!/bin/bash

#This script merges windows of regions that are within certain distances of each other.


#input_fns=$(ls ../data/All_EES_*sorted*.bed)
declare -a distances=("50" "100" "150" "200" "400" "600" "800" "1000" "1500" "2000")
#declare -a chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "ch17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for dist in "${distances[@]}"
do
	echo $dist
	input_files=$(ls '/home/mudrachoudhury/Xiaolab/dsRNA_detection_2020/data/Merged_'$dist'_dist_All_EES_50_window_regions_w_min_3_ES_per_window_using_REDI_2020_'*'_sorted.bed')
	output_file='/home/mudrachoudhury/Xiaolab/dsRNA_detection_2020/data/Merged_'$dist'_dist_All_EES_50_window_regions_w_min_3_ES_per_window_using_REDI_2020_combined.bed'
	if [ -f "$output_file" ] ; then
    		rm "$output_file"
	fi
	cat $input_files >> $output_file
done
