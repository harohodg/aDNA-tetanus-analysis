#!/usr/bin/env bash

#Script to launch bwa/0.7.17 indexing jobs on compute canada infrastructure
#Since they are expected to take less then a minute each we're using parallel instead of sbatch
#Usage <input folder>
#Will search for all .fa files under inputfolder and runs ./bwa_index_job.sh against it

find "$1" -name "*.fa" | parallel --eta -j 4 ./bwa_index_job.sh {}

