#!/usr/bin/env bash

#Script to launch megahit 1.2.9 jobs on Compute Canada Infrastructure
#Usage ./launch_megahit_jobs.sh <input_folder> <output_folder>  <min_mem (gigabytes)> <min_time (minutes)> <num_cpu> <numGPU> <account>
#Replicates the folder structure under input folder in output folder

#If target mem > 80G then entire node will be used

#Has been tested on Beluga
#May need some tweaking for Cedar and Graham

if [ "$#" -ne 7 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./launch_megahit_jobs.sh <input_folder> <output_folder>  <min_mem (gigabytes)> <min_time (minutes)> <num_cpu> <numGPU> <account>"
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


existingJobs=$(sq -h --format="%j" | grep ".*_megahit" | cut -f1 -d'_' | sed 's|-|/|g' | sort)
alreadyRun=$(find "$outputFolder" -name "final.contigs.fa" -printf '%h\n'  | sed "s|$outputFolder||" | sort)
doNotRun=$(sort --unique  <(echo "$existingJobs") <(echo "$alreadyRun") )
#maybeRun=$(find "$inputFolder" -type f -name "*.fastq" -not -empty -printf '%h\n' | sed "s|$inputFolder||" | sort | uniq)
maybeRun=$(find "$inputFolder" -type f -name "*.f*q" -not -empty -printf '%h\n' | sed "s|$inputFolder||" | sort | uniq)

for sourceFolder in $(comm -23 <(echo "$maybeRun") <(echo "$doNotRun") )
do
    inputDirname="$inputFolder"/"$sourceFolder"
    targetFolder="$outputFolder"/"$sourceFolder"
    inputSize=$(du -BM "$inputDirname" | cut -f1 -d'M')
    
    proposedJobMem=$(expr "$inputSize" / 1000)
    jobMEM=$((proposedJobMem>minMem ? proposedJobMem : minMem))G
    
    proposedJobTime=$(expr "$inputSize" / 9 / 60)
    jobTime="0:$((proposedJobTime>minTime ? proposedJobTime : minTime)):0"
    
    inputFiles=$(ls $inputDirname/*.f*q | tr "\n" " " )
    jobName="$(echo $sourceFolder | sed 's|/|-|g')"_megahit
    
    echo sbatch --account="$account" --cpus-per-task="$numCPU" --gres=gpu:"$numGPU" --time="$jobTime" --mem="$jobMEM" --job-name="$jobName" --output=%x-%j.out megahit_job.sh "$inputFiles" "$targetFolder"
done

