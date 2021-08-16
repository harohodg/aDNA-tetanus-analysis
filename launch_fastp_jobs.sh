#!/usr/bin/env bash

#Script to launch a bunch of fastp/0.20.1 jobs with default parameters
#Usage ./launch_fastp_jobs.sh <input_folder> <output_folder> <min_mem (gigabytes)> <min_time (minutes)> <account>
#Finds all folders with .fastq files in them under input folder and recreates the same folde structure under output folder

#Note : fastp 0.20.1 uses a max of 16 threads
#Script currently does not check if jobs are already running

if [ "$#" -ne 5 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./launch_fastp_jobs.sh <input_folder> <output_folder> <min_mem (gigabytes)> <min_time (minutes)> <account>"
    exit 1
fi

inputFolder="$1"
outputFolder="$2"
minMem="$3"
minTime="$4"
account="$5"

jobCPUS=16
jobMEM=$minMem"G"



##Find sraIDs in the input folder that don't have a corresponding *.json file in the output folder
#for sraID in $(comm -23 <(find "$inputFolder" -mindepth 1 -maxdepth 1 -type d -not -empty -exec bash -c 'basename {}' \; | sort) <(find "$outputFolder" -name "*.json" -exec bash -c 'basename $(dirname {})' \; | sort ))
#do
#    inputDirname="$inputFolder"/"$sraID"
#    targetFolder="$outputFolder"/"$sraID"
#    if [ ! -d "$targetFolder" ]; then
#        mkdir -p "$targetFolder"
#    fi
#    
#    inputSize=$(du -cBG "$inputDirname" | tail -n 1 | cut -f1 -d'G')
#    proposedRunTime=$(expr "$inputSize" \* 15 / 60)
#    runTime="0:$((proposedRunTime>minTime ? proposedRunTime : minTime)):0"  
#    echo sbatch --account="$account" --mem="$jobMEM" --ntasks="$jobCPUS" --time="$runTime" --job-name="$sraID""_fastp" --output=%x-%j.out fastp_job.sh "$inputDirname" "$targetFolder"
#done

existingJobs=$(sq -h --format="%j" | grep ".*_$jobNameFormat" | cut -f1 -d'_' | sed 's|-|/|g' | sort)
alreadyRun=$(find "$outputFolder" -type f -name "*.json" | sed "s|$outputFolder||" | sort)
doNotRun=$(sort --unique  <(echo "$existingJobs") <(echo "$alreadyRun") )
maybeRun=$(find "$inputFolder" -type f -name "*.fastq" -not -empty -printf '%h\n' | sed "s|$inputFolder||" | sort | uniq)



for sourceFolder in $( comm -23 <(echo "$maybeRun") <(echo "$doNotRun") )
do
    echo "$sourceFolder"
    inputDirname="$inputFolder"/"$sourceFolder"
    targetFolder="$outputFolder"/"$sourceFolder"
    
    inputSize=$(du -cBG "$inputDirname" | tail -n 1 | cut -f1 -d'G')
    proposedRunTime=$(expr "$inputSize" \* 15 / 60)
    runTime="0:$((proposedRunTime>minTime ? proposedRunTime : minTime)):0" 
    
    jobName="$(echo $sourceFolder | sed 's|/|-|g')"_fastp 
    
    echo sbatch --account="$account" --mem="$jobMEM" --time="$runTime" --job-name="$jobName" --output=%x-%j.out fastp_job.sh "$inputDirname" "$targetFolder"
    
done

