#!/bin/bash

ja=submit_Identify_dsRNAs_step5.job
input_dir="/home/mudrachoudhury/Xiaolab/dsRNA_detection_2020/data"
declare -a distances=("1000")

for n in "${distances[@]}"; do
    input_dir_files=$(ls $input_dir'/Merged_'$n'_EERs/Merged_'$n'_split_chromosomes'*'output.bed')
    for input_file in $input_dir_files; do
        input_file_prefix=$(echo "${input_file%.*}")
        output_file=$input_file_prefix'_w_RNA_fold.bed'
        echo "$input_file $output_file" >> $ja
    done
done
 echo job completed



