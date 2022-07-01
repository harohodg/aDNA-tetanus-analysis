#!/usr/bin/env bash

#Script to relaunch megahit jobs which ran out of time or memory
#$1 = root folder to search
#Assumed to be folder/SRA_ID/SRA_files
#$2 = the results folder for the assembled files
#$3 = the sbatch account to use 
#$4 = the folder with the previous runs .out files
#$5 = the run time for all out of time jobs in minutes
#$6 = the total memory to use for all of memory jobs in MB

memCutoff="80000" #Max memory in MB before we use the whole node
minMem="4096" #Minimum memory to use in MB
minTime="1800" #Seconds


#Need to add max of mem and 4000M

echo "Jobs out of Time"
for fname in $(grep "DUE TO TIME LIMIT" "$4"/*.out | cut -f1 -d':' )
do
    sraID=$(basename "$fname" | cut -f1 -d'_')
    #echo "$sraID"
    
    folder="$1"/"$sraID"
    folder_size=$(du -BM "$folder" | cut -f1 -d'M')
    
    mem_size="$6"
    mem_size=$((mem_size>minMem ? mem_size : minMem))
    if [[ $mem_size -gt $memCutoff ]];then
        cores='32'
        mem='127000M'
        gpu='2'
    else
        cores=16
        mem="$mem_size"M
        gpu='1'
    fi
    

    time=0:"$5":0
    
    echo sbatch --job-name="$sraID"_megahit_job --output=%x-%j.out --account="$3" --ntasks="$cores" --gres=gpu:p100:"$gpu" --time="$time" --mem="$mem" megahit_job.sh $(ls $folder/*.fastq | tr "\n" " " ) "$2"/$(basename $folder)_assembled

done


echo "Jobs out of Memory"
for fname in $(grep "killed by the cgroup out-of-memory" "$4"/*.out | cut -f1 -d':' )
do
    sraID=$(basename "$fname" | cut -f1 -d'_')
    #echo "$sraID"
    
    folder="$1"/"$sraID"
    folder_size=$(du -BM "$folder" | cut -f1 -d'M')
    

    mem_size="$6"
    mem_size=$((mem_size>minMem ? mem_size : minMem))
    if [[ $mem_size -gt $memCutoff ]];then
        cores='32'
        mem='127000M'
        gpu='2'
    else
        cores=16
        mem="$mem_size"M
        gpu='1'
    fi
    

    time=0:"$5":0
    
    echo sbatch --job-name="$sraID"_megahit_job --output=%x-%j.out --account="$3" --ntasks="$cores" --gres=gpu:p100:"$gpu" --time="$time" --mem="$mem" megahit_job.sh $(ls $folder/*.fastq | tr "\n" " " ) "$2"/$(basename $folder)_assembled

done
