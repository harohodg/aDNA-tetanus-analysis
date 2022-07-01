#!/usr/bin/env bash

#Script to run leeHom on Compute Canada infrastructure using GNU parallel
#Usage : parallel_leeHom_job.sh <leeHom source folder> <input_folder> <output folder>


leeHomsrcFolder="$1"
inputFolder="$2"
outputFolder="$3"
outputPrefix="$(basename $outputFolder)_leeHom"

threads=2
command="$leeHomsrcFolder"/leeHom
flags="--verbose -t $threads -fqo $outputFolder/$outputPrefix "

numFiles=$(ls "$inputFolder"/*.fastq | wc -l)
inputFiles=$(ls "$inputFolder"/*.fastq | sort | tr '\n' '|')
inputFile1=$(echo $inputFiles | cut -f1 -d'|')
inputFile2=$(echo $inputFiles | cut -f2 -d'|')

if [[ $numFiles -eq 1 ]]; then
    flags="$flags -fq1 $inputFile1"
#If two input files run PE    
elif [[ $numFiles -eq 2 ]]; then
    flags="$flags --ancientdna -fq1 $inputFile1 -fq2 $inputFile2 "
else
	echo "Incorrect number of fastq files in $inputFolder"
	exit 2
fi


module load StdEnv/2020 gcc/10.2.0

echo "Running : $command $flags"

mkdir -p "$outputFolder"

eval "$command" "$flags" && echo "DONE"
