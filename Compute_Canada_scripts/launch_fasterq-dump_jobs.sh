#!/usr/bin/env bash


#Script to run prefetch SRA datasets through fasterq-dump
#Usage ./launch_fasterq-dump_jobs.sh <input folder> <output folder> <account>
#Assumes only successfully prefetched files are in the input folder

for fname in $(find "$1" -name "*.sra" -exec basename {} \;)
do
    sraID=$(echo "$fname" | cut -f1 -d'.')
    echo "$sraID"
    ./launch_fasterq-dump_job.sh "$1"/"$fname" "$2"/"$sraID" "$3"
done
