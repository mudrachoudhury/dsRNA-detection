#!/bin/bash

# Script Description: 
# This script combines separate chromosome files into a single file for each specified merge distance.

# Define an array of merge distances. Currently, it is set to work with a single distance of 1000.
# Uncomment the next line to process multiple distances: 100, 150, 200, 400, 600, and 1000.
declare -a distances=("1000")
# declare -a distances=("100" "150" "200" "400" "600" "1000")

# Loop through each specified merge distance
for dist in "${distances[@]}"
do
    echo $dist  # Print the current merge distance being processed

    # Construct the path for input files and assign them to a variable
    # This collects all files for a given merge distance from split chromosomes.
    input_files=$(ls '/home/mudrachoudhury/Xiaolab/dsRNA_detection_2020/data/Merged_'$dist'_EERs/Merged_'$dist'_split_chromosomes.chr'*'.output_w_RNA_fold.bed')

    # Define the output file path where the combined results will be stored
    output_file='/home/mudrachoudhury/Xiaolab/dsRNA_detection_2020/data/Merged_'$dist'_EERs/Merged_'$dist'_combined_output_w_RNA_fold.bed'

    # Remove the output file if it already exists to avoid appending to an old file
    rm $output_file

    # Concatenate all input files and redirect the output to the defined output file
    cat $input_files >> $output_file
done
