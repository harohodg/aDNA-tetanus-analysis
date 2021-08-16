#!/usr/bin/env bash

#Script to map sequencing reads using bowtie2/2.4.2 and createing a bam file using samtools/1.12 on compute canada infrastructure
#Usage ./launch_bowtie2_compressed_jobs.sh <input_folder> <output_folder> <target_index_file> <min time (minutes)> <account_to_use>

if [ "$#" -ne 5 ]; then
    echo "Incorrect # of arguments."
    echo "Usage ./launch_bowtie2_compressed_jobs.sh <input_folder> <output_folder> <target_index_file> <min time (minutes)> <account_to_use>"
    exit 1
fi


inputFolder="$1"
outputFolder="$2"
indexSource="$3"
minTime="$4"
account="$5"


module load StdEnv/2020 samtools/1.12

indexName=$(basename "$indexSource")

existingJobs=$(sq -h --format="%j" | grep ".*_bowtie_$indexName" | cut -f1 -d'_' | sort)
alreadyRun=$(find "$outputFolder" -type f -name "*$indexName""_mapped.bam" -exec bash -c 'samtools quickcheck -q {} && basename $(dirname {})' \; | sort)
doNotRun=$(sort --unique  <(echo "$existingJobs") <(echo "$alreadyRun") )
maybeRun=$(find "$inputFolder" -mindepth 1 -maxdepth 1 -type d -exec bash -c 'basename {}' \; | sort)


for sraID in $( comm -23 <(echo "$maybeRun") <(echo "$doNotRun") )
do
    mkdir -p "$outputFolder"/"$sraID"
    
    inputDirname="$inputFolder"/"$sraID"
    targetFolder="$outputFolder"/"$sraID"
    inputSize=$(du -cBG "$inputDirname" | tail -n 1 | cut -f1 -d'G')
    proposedJobTime=$(expr "$inputSize" \* 2)   
    jobTime="0:$((proposedJobTime>minTime ? proposedJobTime : minTime)):0"
    
    echo sbatch --account="$account" --time="$jobTime" --job-name="$sraID""_bowtie_$indexName" --output=%x-%j.out bowtie2_compressed_job.sh "$inputDirname"/*.fastq "$indexSource" "$targetFolder"/"$sraID"_"$indexName"_mapped.bam
done
