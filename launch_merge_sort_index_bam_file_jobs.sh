#!/usr/bin/env bash

#Script for merging, sorting, indexing bam files using samtools/1.12 on compute canada infrastructure
#Usage : ./launch_merge_sort_index_bam_file_jobs.sh <account> <input_folder> <output_folder> <label> <job_time (minutes)> <numCPU> >
#Will merge all *.sorted.bam files under each folder in the input folder and put them in output_folder/foldername_label_merged.bam and output_folder/foldername_label_merged.sorted.bam

#script should be updated to check for running jobs

if [ "$#" -ne 6 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./launch_merge_sort_index_bam_file_jobs.sh <account> <input_folder> <output_folder> <label> <job_time (H:M:S)> <numCPU>"
    exit 1
fi


module load StdEnv/2020 samtools/1.12

account="$1"
inputFolder="$2"
outputFolder="$3"
label="$4"
jobTime="$5"
numCPU="$6"


targetPattern='*.sorted.bam'



for folder in $(find "$inputFolder" -mindepth 1 -maxdepth 1 -type d)
do 
    mergedFile="$outputFolder"/$(basename "$folder")/$(basename "$folder")_"$label"_merged.bam
    sortedFile=$(echo "$mergedFile" | sed 's/.bam/.sorted.bam/')

    jobName="$(basename $folder)"_"$label"_merge_bam
    jobMem="$numCPU"G
    
    #If merged file is not fine then launch merge/sort/index job
    if ! samtools quickcheck -q "$mergedFile"; then
        echo "Launching merge/sort/index job for $folder"   
        inputSize=$( du -cBM $(find "$folder" -name "$targetPattern" ) | tail -n1 | cut -f1 -d'M' )
        jobMem=$(expr "$inputSize" + 30000 )M
             
        echo sbatch --cpus-per-task="$numCPU" --mem="$jobMem" --account="$account" --time="$jobTime" --job-name="$jobName" --output=%x-%j.out ./merge_bam_files.sh "$folder" "$mergedFile" "$targetPattern"
    #If merged file is fine but sorted file is not then launch sort/index job
    elif ! samtools quickcheck -q "$sortedFile"; then
        echo "Launching sort/index job for $folder"
        inputSize=$(du -BM "$mergedFile" | cut -f1 -d'M')
        jobMem=$(expr "$inputSize" + 30000 )M
        
        echo sbatch --cpus-per-task="$numCPU" --mem="$jobMem" --account="$account" --time="$jobTime" --job-name="$jobName" --output=%x-%j.out ./sort_index_bam_file.sh "$mergedFile" "$sortedFile"
    fi
    #Otherwise assume it finished correctly and move on to next folder
    
    
    
#    samtools quickcheck -q "$testFile"
#    if [ "$?" -ne 0 ];then
#        #echo "$targetFile"
#        echo sbatch --cpus-per-task="$numCPU" --mem="$numCPU"G --account="$account" --time="$jobTime" --job-name= --output=%x-%j.out ./merge_bam_files.sh "$folder" "$targetFile" "$targetPattern"
#    fi
done

#find "$inputFolder" -mindepth 1 -maxdepth 1 -type d -exec bash -c "echo {}" \;


