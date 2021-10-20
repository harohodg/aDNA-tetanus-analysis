#!/usr/bin/env bash

#Script for running bwa align samse jobs on compute canada infrastructure
#Usage ./launch_bwa_align_samse_jobs.sh <n> <o> <l> <contigs_folder> <input_folder> <output_folder>  <account>
#Assumes every folder in <input_folder> represents an sraID and that there is a corresponding folder in <contigs_folder> with a bwa index contig.fa file in it 
#Will take every <input_folder>/sraID/fname.fq and run it through the pipeline into <output_folder>/<sra_id>/fname.mappedTo.contig.sorted.bam 

#Should be updated to mirror input folder structure under output folder
#Can still assume that the last folder name is the corresponding contigs folder name

#Revisions:
#October 20, 2021
# - inputFolders checked to make sure they have only one *.fq file in them
# - contigs folder checked to see if there is a corresponding sraID folder
# - contigs/sraID folder checked to see if there is only one *.fa file in them
# Based off code from https://www.ducea.com/2009/03/05/bash-tips-if-e-wildcard-file-check-too-many-arguments/


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

existingJobs=$(sq -h --format="%j" | grep ".*_$jobNameFormat" | cut -f1 -d'_' | sort)
alreadyRun=$(find "$outputFolder" -type f -name "*sorted.bam" -exec bash -c 'samtools quickcheck -q {} && basename $(dirname {})' \; | sort)
doNotRun=$(sort --unique  <(echo "$existingJobs") <(echo "$alreadyRun") )
maybeRun=$(find "$inputFolder" -mindepth 2 -maxdepth 2 -type f -name "*.fq" -not -empty -exec bash -c 'basename $(dirname {})' \; | sort)



for sraID in $( comm -23 <(echo "$maybeRun") <(echo "$doNotRun") )
do
    #echo "$sraID"
    inputFiles=$(ls  "$inputFolder"/"$sraID"/*.fq 2> /dev/null | wc -l)
    contigFiles=$(ls "$contigsFolder"/"$sraID"/*.fa 2> /dev/null | wc -l)
    if [ "$inputFiles" != "1" ]
    then
        echo "Incorrect number of input files found in $inputFolder/$sraID"
    else
        if [ "$contigFiles" != "1" ]
        then
            echo "Incorrect number of contig files found in $contigsFolder/$sraID"
        else
            inputFile=$(ls "$inputFolder"/"$sraID"/*.fq)
            contigs=$(ls "$contigsFolder"/"$sraID"/*.fa)
            outputFile="$outputFolder"/"$sraID"/$(basename "$inputFile" | sed 's/.fq//').mappedTo.$(basename "$contigs" | sed 's/.fa//').sorted.bam
            
            inputSize=$(du -BM "$inputFile" | cut -f1 -d'M')
            jobTime=0:0:"$(expr $inputSize \* 2)"
            
            sbatch --account="$account" --time="$jobTime" --job-name="$sraID"_"$jobNameFormat" --output=%x-%j.out bwa_align_samse_job.sh "$n" "$o" "$l" "$contigs" "$inputFile" "$outputFile"
        fi
    fi
done
