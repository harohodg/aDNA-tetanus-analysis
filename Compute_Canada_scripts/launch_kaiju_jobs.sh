#!/usr/bin/env bash

#Script to launch kaiju/1.7.4 jobs on compute canada infrastructure
#Usage : ./launch_kaiju_jobs.sh <input_folder> <output_folder> <path_to/nodes.dmp> <path_to/kaiju_db_refseq.fmi> <min_mem (GB)> <min_time (minutes)> <account>
#Searches for final_contigs.fa files in input folder and runs kaiju/1.7.4 on them putting the results in same folder structure under output folder named parent_folder_kaiju.out


if [ "$#" -ne 7 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./launch_kaiju_jobs.sh <input_folder> <output_folder> <path_to/nodes.dmp> <path_to/kaiju_db_refseq.fmi> <min_mem (GB)> <min_time (minutes)> <account>"
    exit 1
fi

inputFolder="$1"
outputFolder="$2"
nodesFiles="$3"
refseqDB="$4"
minMem="$5"
minTime="$6"
account="$7"

existingJobs=$(sq -h --format="%j" | grep ".*_kaiju" | cut -f1 -d'_' | sed 's|-|/|g' | sort)
alreadyRun=$(find "$outputFolder" -name "*_kaiju.out" -printf '%h\n'  | sed "s|$outputFolder||" | sort)
doNotRun=$(sort --unique  <(echo "$existingJobs") <(echo "$alreadyRun") )
maybeRun=$(find "$inputFolder" -type f -name "final.contigs.fa" -not -empty -printf '%h\n' | sed "s|$inputFolder||" | sort | uniq)


for sourceFolder in $(comm -23 <(echo "$maybeRun") <(echo "$doNotRun") )
do
    fname="$inputFolder"/"$sourceFolder"/"final.contigs.fa"
    parentFolder=$(basename "$sourceFolder")
    
    targetFile="$outputFolder"/"$sourceFolder"/"$parentFolder"_kaiju.out
    jobName="$parentFolder"_kaiju
    ./launch_kaiju_job.sh "$nodesFiles" "$refseqDB" "$fname" "$targetFile" "$minMem" "$minTime" "$jobName" "$account"
done
