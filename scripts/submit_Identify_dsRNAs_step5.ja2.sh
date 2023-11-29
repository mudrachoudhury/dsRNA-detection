#!/bin/bash
#SBATCH -o Identify_dsRNAs_step5_RNAfold.out
#SBATCH -e Identify_dsRNAs_step5_RNAfold.err
#SBATCH -p long
#SBATCH --mem 32G
#SBATCH -V
#SBATCH -n 4
#SBATCH --exclude=compute-4-5

ja=submit_Identify_dsRNAs_step5_unfinished.job
job_name=Identify_dsRNAs_step5_unfinished

PARMS=($(awk "NR==$SLURM_ARRAY_TASK_ID" $ja))
input_bed_file=${PARMS[0]} #bed file of all EES regions after applying some length cut offs
output_file=${PARMS[1]} #output file of all the RNAfold information and annotations for each EES region


echo -e "Job begin: $(date)\n" 1>>log/$job_name.$SLURM_ARRAY_TASK_ID.out

python Identify_dsRNAs_step5_RNAfold.py --i $input_bed_file --o $output_file 1>>log/$job_name.$SLURM_ARRAY_TASK_ID.out 2>>log/$job_name.$SLURM_ARRAY_TASK_ID.err

echo -e "Job Completed: $(date)\n" 1>>log/$job_name.$SLURM_ARRAY_TASK_ID.out

#sbatch -a 1-6%6 submit_Identify_dsRNAs_step5.ja2.sh
