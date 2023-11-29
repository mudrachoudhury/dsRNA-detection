#!/bin/bash

for file in $(ls '/home/mudrachoudhury/Xiaolab/dsRNA_detection_2020/data/All_EES_50_window_regions_w_min_3_ES_per_window_using_REDI_2020_chr'*'.bed')
do
	echo $file
	sorted_fn_prefix=$(echo $file | cut -f1 -d".")
	sorted_fn=$sorted_fn_prefix'_sorted.bed'	
	~/software/bedtools2/bin/sortBed -i $file > $sorted_fn
done
