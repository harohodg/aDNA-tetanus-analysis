#!/usr/bin/env bash


#SBATCH --mem=10G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=8
#SBATCH --time=0:30:0

#Script to take a align reads using bwa/0.7.17 and then covert to a sorted bam file using 
#samtools/1.12 on compute canada infrastructure
#Usage: sbatch <...> bwa_align_samse_job.sh <n> <o> <l> <prefix> <in.fq> <ouputfile>

if [ "$#" -ne 6 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : bwa_align_samse_job.sh <n> <o> <l> <prefix> <in.fq> <ouputfile>"
    exit 1
fi

n="$1"
o="$2"
l="$3"
prefix="$4"
inputFile="$5"
outputFile="$6"

if [ -n "$SLURM_CPUS_ON_NODE" ]; then
    threads=$SLURM_CPUS_ON_NODE
else
    threads=2
fi


command=$(printf "bwa aln -t %s -n  %s -o %s -l %s %s %s | bwa samse %s /dev/stdin %s | samtools sort --threads %s - -o %s \n" "$threads" "$n" "$o" "$l" "$prefix" "$inputFile" "$prefix" "$inputFile" "$threads" "$outputFile")

echo "Running: $command"

mkdir -p "$(dirname $outputFile)"

module load StdEnv/2020 bwa/0.7.17 samtools/1.12

eval "$command" && echo "Done"
