#!/usr/bin/env bash

#Script for launching seqtk/1.3 jobs on compute canada infrastructure
#Usage <command> <input folder> <output folder> <file suffix> <account>
#Finds all *.fq.gz files under input folder and runs them against seq <command> putting the results in output folder/file_parent/file.file_suffix

command="$1"
inputFolder="$2"
outputFolder="$3"
fileSuffix="$4"
account="$5"

logNameFormat="%x-%j.out"

for fname in $(find "$inputFolder" -name "*.fq.gz" -not -name "*_r[1,2].fq.gz" -not -name "*.fail.fq.gz")
do
    parentFolder="$(basename $(dirname $fname))"
    targetFile="$outputFolder"/"$parentFolder"/$(basename "$fname" | sed "s/fq.gz/$fileSuffix/")
    jobName="$parentFolder"_seqtk
    printf "sbatch --account=%s --job-name=%s --output=%s seqtk_job.sh '%s' %s %s \n" "$account" "$jobName" "$logNameFormat" "$command" "$fname" "$targetFile"
done

