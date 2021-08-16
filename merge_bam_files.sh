#!/usr/bin/env bash
#SBATCH --mem=4G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=4
#SBATCH --time=0:15:0


#Script for merging, sorting, indexing bam files using samtools/1.12 on compute canada infrastructure
#Usage : ./merge_bam_files.sh <input_folder> <target_file> <target_pattern>
#Will merge all bam files matching target pattern in input folder into target file and then sort and index it



inputFolder="$1"
targetFile="$2"
targetPattern="$3"

mkdir -p $(dirname "$targetFile")

./merge_sort_index_bam_files.sh "$targetFile" $(find "$inputFolder" -name "$targetPattern" | tr '\n' ' ')
