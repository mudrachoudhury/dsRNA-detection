#!/bin/bash

# Script Description:
# This script merges windows of genomic regions that are within certain distances of each other.

# Define an array of distances to process. These distances represent the range within which windows will be merged.
declare -a distances=("50" "100" "150" "200" "400" "600" "800" "1000" "1500" "2000")

# Uncomment the following line and edit the script if specific chromosomes are to be processed. Currently, it is commented out.
# declare -a chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "ch17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

# Loop through each specified merge distance
for dist in "${distances[@]}"
do
    echo $dist  # Print the current merge distance being processed

    # Construct the path for input files and assign them to a variable
    # This collects all sorted bed files for a given merge distance.
    input_files=$(ls '/home/mudrachoudhury/Xiaolab/dsRNA_detection_2020/data/Merged_'$dist'_dist_All_EES_50_window_regions_w_min_3_ES_per_window_using_REDI_2020_'*'_sorted.bed')

    # Define the output file path where the combined results for each distance will be stored
    output_file='/home/mudrachoudhury/Xiaolab/dsRNA_detection_2020/data/Merged_'$dist'_dist_All_EES_50_window_regions_w_min_3_ES_per_window_using_REDI_2020_combined.bed'

    # Check if the output file already exists. If it does, remove it to ensure a fresh start.
    if [ -f "$output_file" ] ; then
        rm "$output_file"
    fi

    # Concatenate all input files for the current distance and redirect the output to the defined output file
    cat $input_files >> $output_file
done

