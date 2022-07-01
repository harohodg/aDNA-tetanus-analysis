#!/usr/bin/env bash

#$1 input folder
#$2 output folder
#$3 target index file
#$4 account to use
#$5 the previous jobs folder

minTime=720 #minutes

indexName=$(basename "$3")
for fname in $(grep -L "overall alignment rate" "$5"/*.out)
do
    sraID=$(basename "$fname" | cut -f1 -d'_')
    mkdir -p "$2"/"$sraID"
    inputDirname="$1"/"$sraID"
    inputSize=$(du -cBG "$inputDirname" | tail -n 1 | cut -f1 -d'G')
    runTime=$(expr "$inputSize" \* 10)   
    echo sbatch --account="$4" --time=0:$((runTime>minTime ? runTime : minTime)):0 --job-name="$sraID""_bowtie_$indexName" --output=%x-%j.out bowtie2_compressed_job.sh "$inputDirname"/*.fastq "$3" "$2"/"$sraID"/"$sraID"_"$indexName"_mapped.bam
done
