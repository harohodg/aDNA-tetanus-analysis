#!/usr/bin/env bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=4G
#SBATCH --time=0:30:0


#Script to run fasterq-dump on a prefetched SRA dataset
#$1 the prefetched .sra file
#$2 the target folder to put the results in

module load nixpkgs/16.09 sra-toolkit/2.9.6

if [ -n "$SLURM_TMPDIR" ]; then
    tmpDirectory="$SLURM_TMPDIR"
    cp "$1" "$tmpDirectory"
    sourceFile="$tmpDirectory"/$(basename "$1")
else
    tmpDirectory="./fasterq-dump_tmp_files"
    sourceFile="$1"
fi

fasterq-dump -t "$tmpDirectory" -p --verbose -S "$sourceFile" -O "$2"
