#!/usr/bin/env bash

#Script to launch megahit jobs
#$1 = root folder to search
#Assumed to be folder/SRA_ID/SRA_files
#$2 = the results folder for the assembled files
#$3 = the sbatch account to use 

memCutoff="80000" #Max memory in MB before we use the whole node
minMem="2048" #Minimum memory to use

mkdir -p "$2"

#Need to add max of mem and 4000M

echo "SE Jobs"
#Find and launch all SE jobs
for folder in $(comm -23 <(find $1 -not -name "*1.fastq" -not -name "*2.fastq" -name "*.fastq" -type f -exec dirname {} \; | cut -f3 -d '/' | sort) <(find ../untrimmed_assembled_data/ -type f -name final.contigs.fa -exec dirname {} \; | cut -f3 -d'/' | cut -f1 -d'_' | sort))
do
    folder="$1"/$folder
    
    folder_size=$(du -BM "$folder" | cut -f1 -d'M')
    mem_size=$(expr "$folder_size")
    mem_size=$((mem_size>minMem ? mem_size : minMem))
    if [[ $mem_size -gt $memCutoff ]];then
        cores='32'
        mem='127000M'
        gpu='2'
    else
        cores=16
        mem="$mem_size"M
        gpu='1'
    fi
    
    if [ -n "$4" ]; then
        time="$4"
    else
        time="0:0:$(expr "$folder_size" / 9)"
    fi
    
    echo sbatch --account="$3" --ntasks-per-node="$cores" --gres=gpu:p100:"$gpu" --time="$time" --mem="$mem" megahit_job.sh $(ls $folder/*.fastq | tr "\n" " " ) "$2"/$(basename $folder)_assembled
done

echo "PE Jobs"
#Find and launch all PE jobs
for folder in $(comm -23 <(find $1 -name "*_1.fastq" -type f -exec dirname {} \; | cut -f3 -d '/' | sort) <(find ../untrimmed_assembled_data/ -type f -name final.contigs.fa -exec dirname {} \; | cut -f3 -d'/' | cut -f1 -d'_' | sort))
do
    folder="$1"/$folder

    folder_size=$(du -BM "$folder" | cut -f1 -d'M')
    mem_size=$(expr "$folder_size")
    mem_size=$((mem_size>minMem ? mem_size : minMem))
    if [[ $mem_size -gt $memCutoff ]];then
        cores='32'
        mem='127000M'
        gpu='2'
    else
        cores=16
        mem="$mem_size"M
        gpu='1'
    fi
    
    if [ -n "$4" ]; then
        time="$4"
    else
        time="0:0:$(expr "$folder_size" / 9)"
    fi
    
    sbatch --account="$3" --nodes=1 --ntasks=1 --cpus-per-task=32  --time="$time" --mem=250G megahit_job.sh $(ls $folder/*.fastq | tr "\n" " " ) "$2"/$(basename $folder)_assembled
done
