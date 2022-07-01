#!/usr/bin/env bash

#SBATCH --mem=4G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=2


#Script to run leeHom on Compute Canada infrastructure
#Usage : <leeHom source folder> <output folder> <Output fastq prefix>  <input file1> [<input file2>]


leeHomsrcFolder="$1"
outputFolder="$2"
outputPrefix="$3"

if [ -n "$SLURM_CPUS_ON_NODE" ]; then
    threads=$SLURM_CPUS_ON_NODE
else
    threads=2
fi

command="$leeHomsrcFolder"/leeHom
flags="--verbose -t $threads -fqo $outputFolder/$outputPrefix "

if [[ $# -eq 4 ]]; then
    flags="$flags -fq1 $4"
elif [[ $# -eq 5 ]]; then
    flags="$flags --ancientdna -fq1 $4 -fq2 $5 "   
else
	echo "Incorrect number of parameters"
	echo "Usage : <leeHom source folder> <output folder> <Output fastq prefix>  <input file1> [<input file2>]"
	exit 2
fi



module load StdEnv/2020 gcc/10.2.0

eval "$command" "$flags"
