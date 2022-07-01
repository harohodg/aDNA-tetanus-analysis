#!/usr/bin/env bash

#Script for launching kaiju/1.7.4 job on compute canada infrastructure
#Usage ./launch_kaiju_job.sh <path_to/nodes.dmp> <path_to/kaiju_db_refseq.fmi> <input_contigs> <kaiju_out_file> <min_mem (GB)> <min_time (minutes)> <job name> <account>

nodesFile="$1"
refseqDB="$2"
inputFile="$3"
outputFile="$4"
minMem="$5"
minTime="$6"
jobName="$7"
account="$8"

inputSize=$(du -cBG "$nodesFile" "$refseqDB" "$inputFile" | tail -n1 | cut -f1 -d'G'  )

proposedJobMem=$(expr "$inputSize" + 25)
jobMem=$((proposedJobMem>minMem ? proposedJobMem : minMem))G


proposedRunTime=$(expr "$inputSize" / 1000 / 60  )
runTime="0:$((proposedRunTime>minTime ? proposedRunTime : minTime)):0"



echo sbatch --account="$account" --time="$runTime" --mem="$jobMem" --job-name="$jobName" --output=%x-%j.out kaiju_job.sh "$nodesFile" "$refseqDB" "$inputFile" "$outputFile"
