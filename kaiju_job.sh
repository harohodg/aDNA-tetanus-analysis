#!/usr/bin/env bash

#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks
#SBATCH --cpus-per-task=16

#Script to run kaiju/1.7.4 on Compute Canada infrastructure
#Usage sbatch <...> kaiju_job.sh <path_to/nodes.dmp> <path_to/kaiju_db_refseq.fmi> <input_contigs> <kaiju_out_file>

nodesFiles="$1"
refseqDB="$2"
inputFile="$3"
outputFile="$4"


module load StdEnv/2020  gcc/9.3.0 kaiju/1.7.4

mkdir -p $(dirname "$outputFile")

kaiju -v -t "$nodesFiles" -f "$refseqDB" -i "$inputFile" -o "$outputFile" -z $SLURM_CPUS_ON_NODE && echo "DONE"
