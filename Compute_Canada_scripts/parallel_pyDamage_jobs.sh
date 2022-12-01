#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=parallel_pyDamage_jobs
#SBATCH --output=%x-%j.out

#Version History 
#Version 1.0 : October 19, 2022
#   Functional script with no error checking
#   Currently assumes the files are sitting in scratch

inputRootFolder="$1"
outputRootFolder="$2"
pyDamageApptainerFile="$3"
subFolder="$4" #Set to grouped if you want the data grouped

module load apptainer/1.0

outputRootFolder="${outputRootFolder}/${subFolder:-ungrouped}"; \
pyDamageFlags="--verbose ${subFolder:+--group}" \
outFileLabel=".${inputRootFolder##*/}.${outputRootFolder##*/}_pydamage_results" \
inputRootFolder=$(realpath ${inputRootFolder}  | sed 's|.*scratch|/scratch|') \
outputRootFolder=$(realpath --canonicalize-missing ${outputRootFolder}  | sed 's|.*scratch|/scratch|') \
export inputRootFolder outputRootFolder outFileLabel pyDamageFlags pyDamageApptainerFile 
find "${inputRootFolder}" -name '*.sorted.bam' -printf '%P\n' | parallel \
    'inputFile="${inputRootFolder}/{}";' \
    'outputFolder="${outputRootFolder}/{//}";' \
    'outputFile="${outputFolder}/${outputFolder##*/}${outFileLabel}.csv";' \
    'jobCommand="[[ ! -f $outputFile ]] && apptainer exec -B /home -B /project -B /scratch -B /localscratch $pyDamageApptainerFile pydamage --outdir $outputFolder analyze $pyDamageFlags $inputFile && mv $outputFolder/pydamage_results.csv $outputFile";'\
    'eval $jobCommand' 
unset inputRootFolder outputRootFolder outFileLabel pyDamageFlags pyDamageApptainerFile
