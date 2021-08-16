#!/usr/bin/env bash


#Script to launch a make blast dbs job
#$1 = root folder to search
#Assumed to be folder/SRA_ID_assebled/assembled SRA_files
#$2 = the results folder for the database files
#$3 = time in H:M:S

#Run using bash launch_make_blast_dbs.sh <input_folder> <output_folder>

num_files=$(find "$1" -name final.contigs.fa  | wc -l) 

sbatch --time="$3" make_blast_dbs.sh "$1" "$2"

