#!/bin/bash
#SBATCH -o Identify_dsRNAs_step3.out
#SBATCH -e Identify_dsRNAs_step3.err
#SBATCH -p medium
#SBATCH --mem 12G
#SBATCH -V

ja=Identify_dsRNAs_step3.job
job_name=Identify_dsRNA_step3

#This script merges windows of regions that are within certain distances of each other. Then the script increases the EES regions by 50 bp on each side and merges overlapping regions.
idx=${SLURM_ARRAY_TASK_ID}

PARMS=($(awk "NR==$idx" $ja)) #task ID is the line of the job file (aka job number)
input_fns=${PARMS[0]}

echo -e "Job begin: $(date)\n" 1>>log/$job_name.$idx.out 2>>log/$job_name.$idx.err

#input_fns=$(ls ../data/All_EES_*sorted*.bed)
declare -a distances=("50" "100" "150" "200" "400" "600" "800" "1000" "1500" "2000")

for i in $input_fns; do
    for n in "${distances[@]}"; do
        new_file=$(echo $i | sed -e 's@.*/@@' )
        output_fn="../data/Merged_"$n"_dist_"$new_file
        temp_fn1=$output_fn'_temp1.bed'
        temp_fn2=$output_fn'_temp2.bed'
        echo output file: $output_fn 1>>log/$job_name.$idx.out
        ~/software/bedtools2/bin/mergeBed -d $n -s -i $i > $temp_fn1 2>>log/$job_name.$idx.err #We first merge all regions within a given distance together and overlapping regions
        awk -v OFS='\t' '{print $1, $2-50, $3+50, $4}' $temp_fn1 > $output_fn 2>>log/$job_name.$idx.err #We increase the length of each region by 50bp before and 50bp after the start and end
        #~/software/bedtools2/bin/mergeBed -s -i $temp_fn2 > $output_fn 2>>log/$job_name.$idx.err #We again merge the any regions that are now overlapping
        rm $temp_fn1
        #rm $temp_fn2
        echo Done 1>>log/$job_name.$idx.out
    done
done

echo -e "Job completed: $(date)\n" 1>>log/$job_name.$idx.out 2>>log/$job_name.$idx.err

