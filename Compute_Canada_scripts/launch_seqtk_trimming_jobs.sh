#!/usr/bin/env bash

#SBATCH --mem-per-cpu=1G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=16               # total number of tasks across all nodes
#SBATCH --time=0:30:0


#Script to launch seqtk/1.3 with seq -L? on Compute canada infrastructure
#Usage: ./launch_seqtk_trimming_jobs.sh <input_folder> <output_folder> <trim_length> <search_pattern> <parallels_log_file>

#Searches for "search pattern" files in the input_folder and puts them in outputfolder/..../file.trim_length.fq
#Jobs are run using GNU parallel

#Currently does not check for already running jobs or pre existing jobs

if [ "$#" -ne 5 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./launch_seqtk_trimming_jobs.sh <input_folder> <output_folder> <trim_length> <search_pattern> <parallels_log_file>"
    exit 1
fi


inputFolder="$1"
outputFolder="$2"
trimLength="$3"
searchPattern="$4"
logFile="$5"

if [ -n "$SLURM_CPUS_ON_NODE" ]; then
    threads=$SLURM_CPUS_ON_NODE
else
    threads=2
fi

module load nixpkgs/16.09  gcc/7.3.0 seqtk/1.3

find "$inputFolder" -type f -name "$searchPattern" | parallel --resume-failed --progress -j "$threads" --joblog "$logFile" ./seqtk_trimming_job.sh "$inputFolder" "$outputFolder" "$trimLength" "{}"


