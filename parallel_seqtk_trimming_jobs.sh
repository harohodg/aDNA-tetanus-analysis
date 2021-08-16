#!/usr/bin/env bash

#SBATCH --mem-per-cpu=1G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=16               # total number of tasks across all nodes
#SBATCH --time=0:30:0

#Script to run a collection of short seqtk trimming jobs on Compute Canada infrastructure using GNU parallel
#Usage: sbatch <...> ./parallel_seqtk_trimming_jobs.sh <logfile> '<command1>' '[command 2]' ...

logFile="$1"


if [ -n "$SLURM_CPUS_ON_NODE" ]; then
    threads=$SLURM_CPUS_ON_NODE
else
    threads=2
fi


cat <( for command in "${@:2}";do printf '%s\n' "$command";done) | parallel --dry-run --progress -j "$threads" --joblog "$logFile" {}

