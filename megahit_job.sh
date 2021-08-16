#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16


#Script to run megahit 1.2.9 job on Compute Canada
#Tested on Graham

#Run using 'sbatch [--account=?] [--time=?] [--ntasks=?] [--mem=?] megahit_job.sh input_files results_folder'

module purge
module load StdEnv/2020 megahit/1.2.9


program="megahit -t $SLURM_CPUS_ON_NODE --continue"

#Figure out which version of megahit to run
#If one input file run SE
if [[ $# -eq 2 ]]; then
    input_files="-r $1"
    results_folder="$2"
    
#If two input files run PE    
elif [[ $# -eq 3 ]]; then
    input_files="-1 $1 -2 $2"
    results_folder="$3"
else
	echo "Incorrect number of parameters"
	exit 2
fi

mkdir -p "$(dirname $results_folder)"

#Run megahit
eval "$program $program_switch $input_files -o $results_folder"


