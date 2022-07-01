#!/usr/bin/env bash

#Helper script to run mapDamage jobs using GNU Parallel
#Assumes that input_file looks like <input_folder>/sraID/*.bam
#And that the corresponding reference file is contigs_folder/sraID/*.fa
#Will put the results in output_folder/sraID with the title = sraID

#Usage ./run_mapDamage.sh <singularity_file> <output_folder> <contigs_folder> <input_file>

if [ "$#" -ne 5 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./run_mapDamage.sh <singularity_file> <input_folder> <output_folder> <contigs_folder> <input_file>"
    exit 1
fi


singularityFile="$1"
inputFolder="$2"
outputFolder="$3"
contigsFile="$4"
bamFile="$5"

sourceFolder=$(dirname "$bamFile" | sed "s|$inputFolder||")
label="$sourceFolder"
targetFolder="$outputFolder"/"$sourceFolder"


./mapDamage_job.sh "$singularityFile" "$bamFile" "$contigsFile" "$targetFolder" "$label"
