#!/usr/bin/env bash

#Script to download sra datasets to Compute Canada Server
#Usage download_SRA_datasets.sh <input_file> <output_folder>


module load nixpkgs/16.09 sra-toolkit/2.9.6

inputFile="$1"
outputFolder="$2"

find "$outputFolder" -mindepth 1 -maxdepth 1 -type d -not -empty -exec basename {} \; | sort > already_downloaded.txt
    for sraID in $( comm -23 <(sort "$inputFile") already_downloaded.txt)
    do
        echo "$sraID"
        time fasterq-dump -p -S $sraID -O "$outputFolder"/"$sraID"
    done
rm already_downloaded.txt
