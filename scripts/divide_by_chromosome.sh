#!/bin/bash

dist=1000
input_file='../data/Merged_'$dist'_dist_All_EES_50_window_regions_w_min_3_ES_per_window_using_REDI_2020_combined_300-4000_length_cutoff.bed'
output_dir='../data/Merged_1000_EERs'

mkdir $output_dir

for chr in `cut -f 1 $input_file | sort | uniq`; do
	echo $chr
	grep -w $chr $input_file > $output_dir'/Merged_'$dist'_split_chromosomes.'$chr'.output.bed'
done













