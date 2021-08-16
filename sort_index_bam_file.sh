#!/usr/bin/env bash

#SBATCH --mem=4G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=4
#SBATCH --time=0:15:0


#Script to sort and index bam files
#$1 input file
#$2 sorted file
#Usage sbatch [--account=?] [--time=?] [--cpus-per-task=?] [--mem=?] sort_index_bam_file.sh <input_bam_file> <sorted_bam_file>

#Load required modules
module load StdEnv/2020 samtools/1.12

inputFile="$1"
outputFile="$2"
threadMem="768M"



if [[ $# -ne 2 ]]; then
    echo "Incorrect number of paramabers"
    echo "Usage : sort_index_bam_file.sh <input_bam_file> <sorted_bam_file>"
fi


#echo "Sorting/Indexing $inputFile -> $outputFile" using "$threadMem" per thread


if [ -n "$SLURM_CPUS_ON_NODE" ]; then
    threads="$SLURM_CPUS_ON_NODE"
else
    threads=1
fi

if [ -n "$SLURM_MEM_PER_NODE" ]; then
    totalMem="$SLURM_MEM_PER_NODE"
    threadMem=$(echo "( $totalMem - 250*$threads) / $threads" | bc )M
fi

echo "Deleting $outputFile.tmp.*.bam"
rm "$outputFile".tmp.*.bam

echo "Samtools Sorting $inputFile -> $outputFile using $threads threads and $threadMem memory per thread"
samtools sort  --threads "$threads" -m "$threadMem" -o "$outputFile"  "$inputFile" 

echo "Samtools Indexing $outputFile using $threads threads"
samtools index -@ "$threads" "$outputFile"

echo "DONE"


