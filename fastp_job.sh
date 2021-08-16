#!/usr/bin/env bash

#SBATCH --mem=4G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=16


#Small script to run fastp/0.20.1 on Compute Canada infrastructure
#Usage sbatch .. ... <input_folder> <output_folder> 

inputFolder="$1"
outputFolder="$2"

program="fastp"
trimming_parameters="--thread $SLURM_CPUS_ON_NODE -h $outputFolder/fastp.html -j $outputFolder/fastp.json"

#Check the number of input files
numFiles=$(ls "$inputFolder"/*.fastq | wc -l)
inputFiles=$(ls "$inputFolder"/*.fastq )
if [[ $numFiles -eq 1 ]]; then
    input_files="-i $inputFiles"
    output_files="-o $outputFolder/$(basename "$inputFiles" | sed 's/\(.*\)\..*/\1/')_trimmed.fastq"
#If two input files run PE    
elif [[ $numFiles -eq 2 ]]; then
    file1=$(echo $inputFiles | cut -f1 -d' ')
    file2=$(echo $inputFiles | cut -f2 -d' ')
    input_files="-i $file1 -I $file2"
    results_folder="$3"
    output_files="-o $outputFolder/$(basename $file1 | sed 's/\(.*\)\..*/\1/')_filt.fastq -O $outputFolder/$(basename $file2 | sed 's/\(.*\)\..*/\1/')_filt.fastq"
else
	echo "Incorrect number of fastq files in $inputFolder"
	exit 2
fi



#Load modules
module load StdEnv/2020 fastp/0.20.1

mkdir -p "$outputFolder"

#Now run fastp
eval "$program $trimming_parameters $input_files $output_files"

echo "DONE"


