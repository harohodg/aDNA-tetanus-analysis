#!/usr/bin/env bash

#Script to prefetch SRA datasets in parallel
#$1 the file containing a list of SRA IDs one per line
#$2 the target folder to put the prefetched data in
#$3 if given the number of datasets to download at the same time (default 4)

module load nixpkgs/16.09 sra-toolkit/2.9.6

mkdir -p "$2"

if [ -n "$SLURM_TMPDIR" ]; then
    numJobs="$3"
else
    numJobs="4"
fi

currentTime=$(date +%d-%b-%Y_%I:%M%p)
screen -s "SRA_prefetch_$currentTime" -L -Logfile "$currentTime"_SRA_prefetch.log -d -m bash -c "parallel --eta -j $numJobs ./prefetch_SRA_dataset.sh {} $2 < $1"
