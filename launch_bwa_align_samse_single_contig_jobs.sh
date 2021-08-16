#!/usr/bin/env bash

#Script for running bwa align samse jobs on compute canada infrastructure
#Usage ./launch_bwa_align_samse_jobs.sh <n> <o> <l> <contigs_folder> <input_folder> <output_folder>  <account>
#Finds all *.fq files under input_folder and creates a corresponding folder structure under output folder
 

if [ "$#" -ne 7 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./launch_bwa_align_samse_jobs.sh <n> <o> <l> <contigs_folder> <input_folder> <output_folder> <account>"
    exit 1
fi

module load StdEnv/2020 samtools/1.12

n="$1"
o="$2"
l="$3"
contigsFolder="$4"
inputFolder="$5"
outputFolder="$6"
account="$7"

jobNameFormat="n""$n"o"$o"l"$l"_bwa_align_samse

existingJobs=$(sq -h --format="%j" | grep ".*_$jobNameFormat" | cut -f1 -d'_' | sed 's|-|/|g' | sort)
alreadyRun=$(find "$outputFolder" -type f -name "*sorted.bam" -exec bash -c 'samtools quickcheck -q {} && dirname {} ' \; | sed "s|$outputFolder||" | sort)
doNotRun=$(sort --unique  <(echo "$existingJobs") <(echo "$alreadyRun") )
maybeRun=$(find "$inputFolder" -type f -name "*.fq" -not -empty -printf '%h\n' | sed "s|$inputFolder||" | sort)


for sourceFolder in $( comm -23 <(echo "$maybeRun") <(echo "$doNotRun") )
do
    echo "$sourceFolder"
    inputFile=$(ls "$inputFolder"/"$sourceFolder"/*.fq)
    contigs=$(ls "$contigsFolder"/*.fa)
    outputFile="$outputFolder"/"$sourceFolder"/$(basename "$inputFile" | sed 's/.fq//').mappedTo.$(basename "$contigs" | sed 's/.fa//').sorted.bam
    
    inputSize=$(du -BM "$inputFile" | cut -f1 -d'M')
    jobTime=0:0:"$(expr $inputSize \* 4)"
    jobName="$(echo $sourceFolder | sed 's|/|-|g')"_"$jobNameFormat"
    
    echo sbatch --account="$account" --time="$jobTime" --job-name="$jobName" --output=%x-%j.out bwa_align_samse_job.sh "$n" "$o" "$l" "$contigs" "$inputFile" "$outputFile"
    
done

