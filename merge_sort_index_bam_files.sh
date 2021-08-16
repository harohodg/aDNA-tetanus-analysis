#!/usr/bin/env bash

#Script to take n bam files and merge them into one using samtools/1.12 on compute canada infrastructure

echo "Merging bam files : ${@:2} -> $1"

module load StdEnv/2020 samtools/1.12

if [ -n "$SLURM_CPUS_ON_NODE" ]; then
    threads="$SLURM_CPUS_ON_NODE"
else
    threads=1
fi

samtools merge -f --threads "$threads" "$1" "${@:2}"
./sort_index_bam_file.sh "$1" $(echo "$1" | sed 's/.bam/.sorted.bam/')
