#!/bin/bash
#SBATCH -o Identify_dsRNAs_step2.out
#SBATCH -e Identify_dsRNAs_step2.err
#SBATCH -p long
#SBATCH --mem 50G
#SBATCH -V


ja=submit_identify_dsRNAs_step2_chromsomes_unfinished_updated.job
log_file=Identify_dsRNAs_step2_updated
idx=${SLURM_ARRAY_TASK_ID}

PARMS=($(awk "NR==$idx" $ja))
chromosome=${PARMS[0]}
input_file=../data/REDI_portal_sites_2020_$chromosome.txt
window_size=50
minimum_number_editing_sites=3

echo -e "Job begin: $(date)\n" 1>>log/$log_file.$chromosome.out 2>>log/$log_file.$chromosome.err

python Identify_dsRNAs_step2_writeChunks_updated.py -i $input_file -s $window_size -n $minimum_number_editing_sites -c $chromosome 1>>log/$log_file.$chromosome.out 2>>log/$log_file.$chromosome.err

echo -e "Job completed: $(date)\n" 1>>log/$log_file.$chromosome.out 2>>log/$log_file.$chromosome.err
