#!/bin/bash
#SBATCH --mem=4gb
#SBATCH --nodes=1
#SBATCH --ntasks=2

#Script to run trimmomatic 0.39 job on Compute Canada
#Tested on Graham

#Run using 'sbatch [--account=?] [--time=?] [--ntasks=?] [--mem=?] trimmomatic_job.sh input_files results_folder'

module purge
module load StdEnv/2020 trimmomatic/0.39


program="java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar"
trimming_parameters="AVGQUAL:20 MINLEN:50 CROP:50"

#Figure out which version of trimmomatic to run
#If one input file run SE
if [[ $# -eq 2 ]]; then
    program_switch="SE -threads $SLURM_CPUS_ON_NODE"
    input_files="$1"
    results_folder="$2"
    output_files="$results_folder/$(basename "$1" | sed 's/\(.*\)\..*/\1/')_trimmed.fastq"
    
#If two input files run PE    
elif [[ $# -eq 3 ]]; then
    program_switch="PE -threads $SLURM_CPUS_ON_NODE"
    input_files="$1 $2"
    results_folder="$3"
    output_files="$results_folder/$(basename $1 | sed 's/\(.*\)\..*/\1/')_filt.fastq $results_folder/$(basename $1 | sed 's/\(.*\)\..*/\1/')_unpair.fastq $results_folder/$(basename $2 | sed 's/\(.*\)\..*/\1/')_filt.fastq $results_folder/$(basename $2 | sed 's/\(.*\)\..*/\1/')_unpair.fastq"
else
	echo "Incorrect number of parameters"
	exit 2
fi


#Create results folder
mkdir -p $results_folder


#Run trimmomatic
eval "$program $program_switch $input_files $output_files $trimming_parameters"


