#!/usr/bin/env bash

#Script to launch megahit 1.2.9 jobs on Compute Canada Infrastructure
#Usage ./launch_megahit_coassembly_jobs.sh <input_folder> <output_folder>  <min_mem (gigabytes)> <min_time (minutes)> <num_cpu> <numGPU> <account>
#Will search for all folders in input_folder that have only 1 subfolder
#Results go in output_folder/folder/folder_coassembly


#Has been tested on Beluga
#Will need some tweaking for Cedar and Graham

if [ "$#" -ne 7 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./launch_megahit_coassembly_jobs.sh <input_folder> <output_folder>  <min_mem (gigabytes)> <min_time (minutes)> <num_cpu> <numGPU> <account>"
    exit 1
fi


maxMem="100"

inputFolder="$1"
outputFolder="$2"
minMem="$3"
minTime="$4"
numCPU="$5"
numGPU="$6"
account="$7"

#Figure out which biosamples are already running
existingJobs=$(sq -h --format="%j" | grep ".*_megahit_coassembly" | cut -f1 -d'_' | sort)
alreadyRun=$(find "$outputFolder" -name "final.contigs.fa" -exec bash -c 'basename $(dirname {}) | cut -f1 -d"_"' \; | sort)
doNotRun=$(sort --unique  <(echo "$existingJobs") <(echo "$alreadyRun") )
maybeRun=$(uniq -d <(find "$inputFolder" -mindepth 2 -maxdepth 2 -type d -exec bash -c 'dirname {} | sort ' \; | sort ) | rev | cut -f1 -d'/' | rev )


for bioSample in $(comm -23 <(echo "$maybeRun") <(echo "$doNotRun") )
do
    echo "$bioSample"
    jobFiles=$(find "$inputFolder"/"$bioSample" -name "*.fastq" )
    inputSize="0"
    for fname in $(echo "$jobFiles")
    do
        inputSize=$(expr "$inputSize" + $(du -BG "$fname" | cut -f1 -d'G'))
    done
    proposedJobMem=$(expr "$inputSize")
    jobMem=$((proposedJobMem>minMem ? proposedJobMem : minMem))G
    
    proposedJobTime=$(expr "$inputSize" \* 1000 / 9 / 60)
    jobTime="0:$((proposedJobTime>minTime ? proposedJobTime : minTime)):0"
    
    jobFiles=$(echo "$jobFiles" | sort | tr '\n' ' ')
    
    echo sbatch  --account="$account" --ntasks-per-node="$numCPU" --gres=gpu:"$numGPU" --time="$jobTime" --mem="$jobMem" --job-name="$bioSample"_megahit_coassembly --output=%x-%j.out megahit_coassembly.sh "$outputFolder"/"$bioSample"/"$bioSample"_coassembly $jobFiles

done

