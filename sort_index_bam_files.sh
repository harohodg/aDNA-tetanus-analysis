#!/usr/bin/env bash

#Script to sort and index bam files using samtools/1.12 on compute canada infrastructure
#Usage ./launch_sort_index_bam_files.sh <input_folder> <account_to_use>

module load StdEnv/2020 samtools/1.12


inputFolder="$1"
account="$2"

runTime=30 #minutes

existingJobs=$(sq -h --format="%j" | grep ".*_samtools_sort_index" | cut -f1 -d'_' | sort)
alreadyRun=$(find "$inputFolder" -type f -name "*sorted.bam" -exec bash -c "samtools quickcheck -q {} && echo {} | sed 's/sorted.bam/bam/'" \; | sort)
doNotRun=$(sort --unique  <(echo "$existingJobs") <(echo "$alreadyRun") )
maybeRun=$(find "$inputFolder" -name "*.bam" -not -name "*.sorted.bam" | sort)



for bamFile in $( comm -23 <(echo "$maybeRun") <(echo "$doNotRun") )
do
    #echo "$bamFile"
    sraID=$(basename "$bamFile" | cut -f1 -d'_')
    targetFile=$(echo "$bamFile" | sed 's/.bam/.sorted.bam/' )
    echo sbatch --account="$account" --time=0:"$runTime":0 --job-name="$sraID"_samtools_sort_index --output=%x-%j.out sort_index_bam_file.sh "$bamFile" "$targetFile"
done
