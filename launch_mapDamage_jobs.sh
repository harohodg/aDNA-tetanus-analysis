#!/usr/bin/env bash

#SBATCH --mem-per-cpu=500M
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=16               # total number of tasks across all nodes
#SBATCH --time=0:30:0


#Script for running many smaller mapDamage jobs Compute canada infrastructure
#Usage: ./aunch_mapDamage_jobs.sh <singularity_file> <input_folder> <output_folder> <contigs_folder> <account> <job time (minutes)>

#Searches for all input_folder/sraID/*sorted.bam files that don't have a corresponding output_folder/sraID/*.pdf



if [ "$#" -ne 6 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./mapDamage_jobs.sh <singularity_file> <input_folder> <output_folder> <contigs_folder> <account> <job time (minutes)>"
    exit 1
fi

module load StdEnv/2020 samtools/1.12


singularityFile="$1"
inputFolder="$2"
outputFolder="$3"
contigsFolder="$4"
account="$5"
jobTime=0:"$6":0

if [ -n "$SLURM_CPUS_ON_NODE" ]; then
    threads=$SLURM_CPUS_ON_NODE
else
    threads=2
fi

alreadyRun=$(find "$outputFolder" -mindepth 2 -maxdepth 2 -name "*.pdf" -exec bash -c 'basename $(dirname {})' \; | sort | uniq)
maybeRun=$(find "$inputFolder" -mindepth 2 -maxdepth 2 -name "*sorted.bam" -not -name "*bam.tmp.*.bam" -exec bash -c 'samtools quickcheck -q {} && basename $(dirname {})' \; | sort | uniq)

for sraID in $(comm -23 <(echo "$maybeRun") <(echo "$alreadyRun") )
do
    bamFile=$(ls "$inputFolder"/"$sraID"/*.sorted.bam)
    referenceContigs=$(ls "$contigsFolder"/"$sraID"/*.fa)
    echo sbatch --account="$account" --time="$jobTime" --job-name="$sraID"_mapDamage --output=%x-%j.out mapDamage_job.sh "$singularityFile" "$bamFile" "$referenceContigs" "$outputFolder"/"$sraID" $sraID
done

