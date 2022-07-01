#!/usr/bin/env bash

#SBATCH --mem-per-cpu=500M
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=16       
#SBATCH --time=0:30:0


#Script for running many smaller mapDamage jobs Compute canada infrastructure
#Usage: sbatch <....> mapDamage_jobs.sh <singularity_file> <input_folder> <output_folder> <contigs_folder>

#Searches for all *.sorted.bam files under input folder and the checks for a corresponding output_folder/..../*.pdf 
#Will replicate folder structure under input_folder under output_folder
#If given a contigs file uses it with all jobs
#If given a contigs folder assumes a corresponding layout to input_folder with the corresponding contigs_file in the same 
#place as each input file

if [ "$#" -ne 4 ]; then
    echo "Incorrect # of arguments."
    echo "Usage : ./mapDamage_jobs.sh <singularity_file> <input_folder> <output_folder> <contigs_folder/contigs_file> "
    exit 1
fi

module load StdEnv/2020 samtools/1.12


singularityFile="$1"
inputFolder="$2"
outputFolder="$3"
contigs="$4"


if [ -n "$SLURM_CPUS_ON_NODE" ]; then
    threads=$SLURM_CPUS_ON_NODE
else
    threads=2
fi

#alreadyRun=$(find "$outputFolder" -mindepth 2 -maxdepth 2 -name "*.pdf" -exec bash -c 'basename $(dirname {})' \; | sort | uniq)
#maybeRun=$(find "$inputFolder" -mindepth 2 -maxdepth 2 -name "*.bam" -not -name "*.tmp.*.bam" -exec bash -c 'samtools quickcheck -q {} && basename $(dirname {})' \; | sort | uniq)

alreadyRun=$(find "$outputFolder" -name "*.pdf" -printf '%h\n'  | sed "s|$outputFolder||" | sort)
maybeRun=$(find "$inputFolder" -type f -name "*sorted.bam" -not -name "*.bam.tmp.*.bam" -exec bash -c 'samtools quickcheck -q {} && dirname {} ' \; | sed "s|$inputFolder||" | sort)


if [[ -d "$contigs" ]];then
    comm -23 <(echo "$maybeRun") <(echo "$alreadyRun") | parallel --eta -j "$threads" bash -c '"'./mapDamage_job.sh $singularityFile '$( ls '$inputFolder'/{}/*.sorted.bam)' '$(ls '$contigs'/{}/*.fa)' $outputFolder/{}/ '$(basename {})' '"'
else #Otherwise we assume it's a file
    comm -23 <(echo "$maybeRun") <(echo "$alreadyRun") | parallel --eta -j "$threads" bash -c '"'./mapDamage_job.sh $singularityFile '$( ls '$inputFolder'/{}/*.sorted.bam)' $contigs $outputFolder/{}/ '$(basename {})' '"'
fi

#mapDamage_job.sh <sif_file> <input_bam_file> <reference_fa_file> <results_folder> <plots title>"

#comm -23 <(echo "$maybeRun") <(echo "$alreadyRun") | xargs -I {} bash -c "ls $inputFolder/{}/*.sorted.bam" | parallel --dry-run -j "$threads" bash -c '"'./run_mapDamage.sh $singularityFile $inputFolder $outputFolder $contigsFile {}'"' 


#printf 'hello %s\n' "'hello'"

#mapDamage_job.sh ~/scratch/tetanus_aDNA_analysis/singularity_containers/mapDamage.sif $(printf "'--merge-reference-sequences --no-stats'")  "$bamFile" "$referenceContigs" "$folder"/"$sraID" "$sraID"
