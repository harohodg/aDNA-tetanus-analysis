#!/usr/bin/env bash


#Script to run seqtk/1.3 with seq -L? on Compute canada infrastructure
#Usage: ./seqtk_trimming_job.sh <input_folder> <output_folder> <trim_length> <inputfile>
#Will take 

inputFolder="$1"
outputFolder="$2"
trimLength="$3"
inputFile="$4"
outputFile=$(echo "$inputFile" | sed 's|'"$inputFolder"'|'"$outputFolder"'|' | sed s/fq.gz/"$trimLength".fq/)

mkdir -p $(dirname "$outputFile")
eval $(printf "./seqtk_job.sh 'seq -L%s' %s %s" "$trimLength" "$inputFile" "$outputFile") 
