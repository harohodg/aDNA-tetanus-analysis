#!/usr/bin/env bash

#Script to launch leeHom jobs on Compute Canada infrastructure
#Usage : ./launch_leeHom_jobs.sh <leeHom source folder> <input folder> <output folder> <min_time (minutes)> <account>
#Runs leeHom against every fastq file under input folder and puts the results in output folder/fastq parent folder

if [ "$#" -ne 5 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./launch_leeHom_jobs.sh  <leeHom source folder> <input folder> <output folder> <min_time (minutes)> <account>"
    exit 1
fi

srcFolder="$1"
inputFolder="$2"
outputFolder="$3"
minTime="$4"
account="$5"

alreadyRun=$(find "$outputFolder" -type f -name "*.fq.gz" -exec bash -c 'basename $(dirname {})' \; | sort | uniq)
maybeRun=$(find "$inputFolder" -type f -name "*.fastq" -exec bash -c 'basename $(dirname {})' \; | sort | uniq)


for parentFolder in $( comm -23 <(echo "$maybeRun") <(echo "$alreadyRun") )
do
    folder=$(find "$inputFolder" -name "$parentFolder" -type d)
    echo "$folder"
#    #parentFolder=$(basename "$folder")
    jobName="$parentFolder"_leeHom
    resultsFolder="$outputFolder"/"$parentFolder"
    
    
    inputSize=$(du -cBM "$folder"/*.fastq | tail -n1 | cut -f1 -d'M')
    jobMem="1G"
    proposedJobTime=$(expr $inputSize / 2 / 60 )
    jobTime="0:$((proposedJobTime>minTime ? proposedJobTime : minTime)):0"
        
    mkdir -p "$resultsFolder"
    echo sbatch --time="$jobTime" --mem="$jobMem" --account="$account" --job-name="$jobName" --output=%x-%j.out leeHom_job.sh "$srcFolder" "$resultsFolder" "$jobName" "$folder"/*.fastq

done

