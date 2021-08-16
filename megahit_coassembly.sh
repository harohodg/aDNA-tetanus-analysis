#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16


#Script to run megahit coassembly 1.2.9 job on Compute Canada
#Does not check to make sure there are the correct number of inputs

#Run using 'sbatch [--account=?] [--time=?] [--ntasks=?] [--mem=?] megahit_coassembly.sh <results_folder> <input_files>'

module load StdEnv/2020 megahit/1.2.9

if [ -n "$SLURM_CPUS_ON_NODE" ]; then
    threads="$SLURM_CPUS_ON_NODE"
else
    threads=16
fi

results_folder="$1"

program="megahit -t $threads --continue"

#Figure out which version of megahit to run
#If files end with _1.gz then assume PE otherwise SE


for fname in "${@:2}"
do
    if [[ "$fname" == *_1_filt.fastq ]];then
        input1="$input1""$fname",
    elif [[ "$fname" == *_2_filt.fastq ]];then
        input2="$input2""$fname",
    elif [[ "$fname" == *_trimmed.fastq ]];then 
        input="$input""$fname",
    else
        echo "$fname does not seem to be SE or PE"
        exit 1
    fi
done

if [ -n "$input" ]; then
    input_files="-r ${input::-1}"
fi
if [ -n "$input1" ] && [ -n "$input2" ];then
    input_files="$input_files -1 ${input1::-1} -2 ${input2::-1}"
fi

#Run megahit
eval "$program  $input_files -o $results_folder"


