#!/usr/bin/env bash

#SBATCH --mem=1G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --time=0:5:0


#Script for running bwa/0.7.17 on compute canada infrastructure
#Usage <input file>

module load StdEnv/2020 bwa/0.7.17

bwa index "$1" && echo "Done"

