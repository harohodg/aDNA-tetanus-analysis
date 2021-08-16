#!/usr/bin/env bash

#Script to launch trimmomatic jobs
#$1 = root folder to search
#Assumed to be folder/SRA_ID/SRA_files
#$2 = the results folder for the trimmed files
#3 = the sbatch account to use 



min_SE_time=5  #minutes
max_SE_time=20 #minutes
max_SE_GB=100

min_PE_time=15  #minutes
max_PE_time=200 #minutes
max_PE_GB=300

echo "SE Jobs"
#Find and launch all SE jobs
for folder in $(find $1 -not -name "*_1.fastq" -not -name "*_2.fastq" -name "*.fastq" -type f -exec dirname {} \;)
do
    folder_size=$(du -BG "$folder" | cut -f1 -d'G')
    #echo "$folder"
    sbatch --account=$3 --time=00:$(./get_job_time.sh $min_SE_time $max_SE_time $max_SE_GB $folder_size):00 trimmomatic_job.sh $(ls $folder/*.fastq | tr "\n" " ") $2$(basename $folder)
done

echo "PE Jobs"
#Find and launch all PE jobs
for folder in $(find $1 -name "*_1.fastq" -type f -exec dirname {} \;)
do
    folder_size=$(du -BG "$folder" | cut -f1 -d'G')
    #echo "$folder"
    sbatch --account=$3 --time=00:$(./get_job_time.sh $min_PE_time $max_PE_time $max_PE_GB $folder_size):00 trimmomatic_job.sh $(ls $folder/*.fastq | tr "\n" " ") $2$(basename $folder)
done
