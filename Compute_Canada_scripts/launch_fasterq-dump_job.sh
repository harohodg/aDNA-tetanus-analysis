#!/usr/bin/env bash

#Script to sbatch launch a fasterq-dump job
#Usage ./launch_fasterq-dump_job.sh <input prefetched sra data> <output folder> <account>

minRunTime=5 #Minutes

inputSize=$(du -BG "$1" | cut -f1 -d'G')
scratchSize=$(expr "$inputSize" \* 1000 \* 6)
minRunTime
echo sbatch --account="$3" --tmp="$scratchSize" --time=0:"$runTime":0 --job-name=$(basename "$1" | cut -f1 -d'.') --output=%x-%j.out fasterq-dump_job.sh "$1" "$2"
