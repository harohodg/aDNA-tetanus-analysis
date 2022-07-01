#!/usr/bin/env bash

#SBATCH --mem=1G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --time=0:5:0


#Script for running mapDamage in a Singularity container on Compute Canada infrastructure
#Usage: sbatch <...> mapDamage_job.sh <sif_file>  <input_bam_file> <reference_fa_file> <results_folder> <plots title>

if [ "$#" -ne 5 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : mapDamage_job.sh <sif_file> <input_bam_file> <reference_fa_file> <results_folder> <plots title>"
    exit 1
fi

sifFile="$1"
mapDamageParameters="--merge-reference-sequences --no-stats"
inputFile="$2"
referenceFile="$3"
resultsFolder="$4"
plotsTitle="$5"

resultsFolderParent=$(dirname "$resultsFolder")
folderName=$(basename "$resultsFolder")
inputFile=$(realpath $inputFile | sed 's|/lustre04||')

mkdir -p "$resultsFolderParent"

module load singularity/3.7


command=$(printf "cd /data; mapDamage %s -i '%s' -r '%s' --title='%s' --folder='%s' \n" "$mapDamageParameters"  "$inputFile"  "$referenceFile" "$plotsTitle" "$folderName"  )
echo "Running: $command"

singularity run -B /home -B /project -B /scratch -B /localscratch -B "$resultsFolderParent":/data "$sifFile" bash -c "$command"

#singularity shell -B /home -B /project -B /scratch -B /localscratch -B "$resultsFolder":/data "$sifFile"
