#makes sure all step5 jobs completed

original_job=submit_Identify_dsRNAs_step5.job
ja=submit_Identify_dsRNAs_step5_unfinished.job
##############DONE CHANGING PARAMETERS###############

rm $ja #remove old unfinished job file
lines=$(awk 'BEGIN{number=0} {number += 1} END{print number}' $original_job)
#echo $lines
for ((i=1;i<=$lines; i+=1))
do
    log_file=log/Identify_dsRNAs_step5.${i}.out
    #err_file=$log_file.${i}.err
    if ! [ -f "$log_file" ]
    then
        echo $log_file does not exist
        awk "NR==$i {print \$0 \" $i\"}" $original_job >> $ja
        continue
    elif ! grep "Step 5 job completed" $log_file >/dev/null
    then
        echo $log_file job was not completed
        awk "NR==$i {print \$0 \" $i\"}" $original_job >>$ja
        continue
    fi
done
echo job completed checking
