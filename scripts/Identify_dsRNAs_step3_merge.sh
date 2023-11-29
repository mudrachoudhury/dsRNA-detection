#!/bin/bash
#SBATCH -o Identify_dsRNAs_step3.out # Standard output will be written to Identify_dsRNAs_step3.out
#SBATCH -e Identify_dsRNAs_step3.err # Standard error will be written to Identify_dsRNAs_step3.err
#SBATCH -p medium                # Specifies the partition (or queue) to submit the job
#SBATCH --mem 12G                # Memory requirement: allocate 12G of RAM for this job
#SBATCH -V                       # Exports the current environment to the job environment

# Variable definitions
ja=Identify_dsRNAs_step3.job     # Job array file name
job_name=Identify_dsRNA_step3    # Name of the job

# Script description:
# This script merges windows of regions that are within certain distances of each other. 
# Then the script increases the EES regions by 50 bp on each side and merges overlapping regions.

# Extracting the task ID from the SLURM array job
idx=${SLURM_ARRAY_TASK_ID}

# Reading parameters from the job array file
PARMS=($(awk "NR==$idx" $ja))    # Extract the line corresponding to the task ID from the job file
input_fns=${PARMS[0]}            # Assign the first parameter from the line as input filenames

# Logging start time and job information
echo -e "Job begin: $(date)\n" 1>>log/$job_name.$idx.out 2>>log/$job_name.$idx.err

# Declare an array of distances for merging windows
declare -a distances=("50" "100" "150" "200" "400" "600" "800" "1000" "1500" "2000")

# Loop over each input file and perform operations
for i in $input_fns; do
    for n in "${distances[@]}"; do
        # Preparing filenames for output and temporary files
        new_file=$(echo $i | sed -e 's@.*/@@' )
        output_fn="../data/Merged_"$n"_dist_"$new_file
        temp_fn1=$output_fn'_temp1.bed'
        
        # Log the output file name
        echo output file: $output_fn 1>>log/$job_name.$idx.out

        # Merge regions within a given distance using bedtools
        ~/software/bedtools2/bin/mergeBed -d $n -s -i $i > $temp_fn1 2>>log/$job_name.$idx.err

        # Increase the length of each region by 50bp at the start and end
        awk -v OFS='\t' '{print $1, $2-50, $3+50, $4}' $temp_fn1 > $output_fn 2>>log/$job_name.$idx.err

        # Clean up temporary files
        rm $temp_fn1

        # Log completion of processing for the current distance
        echo Done 1>>log/$job_name.$idx.out
    done
done

# Log job completion time
echo -e "Job completed: $(date)\n" 1>>log/$job_name.$idx.out 2>>log/$job_name.$idx.err
