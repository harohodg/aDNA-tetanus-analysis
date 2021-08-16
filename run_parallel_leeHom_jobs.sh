#!/usr/bin/env bash

#SBATCH --mem=16G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=16


#Script to run leeHom jobs on Compute Canada infrastructure using parallel
#Usage : sbatch <....> launch_leeHom_jobs.sh <leeHom source folder> <input folder> <output folder> 
#Runs leeHom against every fastq file under input folder and puts the results in output folder with same folder structure

if [ "$#" -ne 3 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./launch_leeHom_jobs.sh  <leeHom source folder> <input folder> <output folder>"
    exit 1
fi

srcFolder="$1"
inputFolder="$2"
outputFolder="$3"

alreadyRun=$(find "$outputFolder" -name "*_leeHom.fq.gz" -printf '%h\n'  | sed "s|$outputFolder||" | sort)
doNotRun=$(sort --unique  <(echo "$existingJobs") <(echo "$alreadyRun") )
maybeRun=$(find "$inputFolder" -type f -name "*.fastq" -not -empty -printf '%h\n' | sed "s|$inputFolder||" | sort | uniq)


#for parentFolder in $( comm -23 <(echo "$maybeRun") <(echo "$alreadyRun") )
#do
#    folder=$(find "$inputFolder" -name "$parentFolder" -type d)
#    echo "$folder"
##    #parentFolder=$(basename "$folder")
#    jobName="$parentFolder"_leeHom
#    resultsFolder="$outputFolder"/"$parentFolder"
#    
#    
#    inputSize=$(du -cBM "$folder"/*.fastq | tail -n1 | cut -f1 -d'M')
#    jobMem="1G"
#    proposedJobTime=$(expr $inputSize / 2 / 60 )
#    jobTime="0:$((proposedJobTime>minTime ? proposedJobTime : minTime)):0"
#        
#    mkdir -p "$resultsFolder"
#    echo sbatch --time="$jobTime" --mem="$jobMem" --account="$account" --job-name="$jobName" --output=%x-%j.out leeHom_job.sh "$srcFolder" "$resultsFolder" "$jobName" "$folder"/*.fastq

#done
comm -23 <(echo "$maybeRun") <(echo "$alreadyRun") | parallel --will-cite --eta ./parallel_leeHom_job.sh "$srcFolder" "$inputFolder"/{} "$outputFolder"/{}
