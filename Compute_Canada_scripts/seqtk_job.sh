#!/usr/bin/env bash

#SBATCH --mem-per-cpu=1G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --time=0:5:0

#Script for running seqtk/1.3 on compute canada infrastructure
#Usage <command> <input file> <output file>


module load nixpkgs/16.09  gcc/7.3.0 seqtk/1.3

printf "Running seqtk %s %s > %s \n" "$1" "$2" "$3"


mkdir -p $(dirname "$3")
eval "seqtk" "$1" "$2" > "$3"
echo "DONE"
