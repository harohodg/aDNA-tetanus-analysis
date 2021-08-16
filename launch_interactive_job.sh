#!/usr/bin/env bash

#Script for launching an interactive job on Compute Canada infrastructure with a log file
#Usage : ./launch_interactive_job.sh <job_name> <log_file> <num_cpu> <total mem> <run time>
#Eg    : ./launch_interactive_job.sh coverage_plots coverage_plots_$(date +"%B_%d_%Y").log 1 12G 1:30:0

if [ "$#" -ne 5 ]; then
    echo 'Incorrect usage.'
    echo 'Usage   : ./launch_interactive_job.sh <job_name> <log_file> <num_cpu> <total mem> <run time>'
    echo 'Example : ./launch_interactive_job.sh coverage_plots coverage_plots_$(date +"%B_%d_%Y").log 1 12G 1:30:0'
    exit 2
fi


jobName="$1"
logFile="$2"
numCPU="$3"
totalMem="$4"
runTime="$5"


screen -S "$jobName" -L -Logfile "$logFile"  salloc --account="def-acdoxey" --nodes=1 --tasks=1 --cpus-per-task="$numCPU" --mem="$totalMem" --time="$runTime" --job-name="$jobName"

