#!/usr/bin/env bash

#SBATCH --mem=2G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=16
#SBATCH --time=0:15:0


#Script to map reads to a reference genome using bowtie2
#If three parameters, assumed to be unpaired reads
#$1 = input file
#$2 = index
#$3 = output BAM file
#If four parameters, assumed to be paired reads
#$1 = input file #1
#$2 = input file #2
#$3 = index
#$4 = output BAM file
#Usage sbatch [--account=?] [--time=?] [--cpus-per-task=?] [--mem=?] bowtie2_job.sh input_file1.fna <input_file2.fna> index SAM file


#Note : With the current defaults phred33 scoring is assumed
#Should possibly be updated to sort and index after bam file is created

#Load required modules
module load StdEnv/2020 bowtie2/2.4.2 samtools/1.12




program="bowtie2"
program_parameters="-p 20 --local -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --threads $SLURM_CPUS_ON_NODE"

if [[ $# -eq 3 ]]; then
    bowtie_input="-x $2 -U $1"
    samtools_output="-o $3"  
elif [[ $# -eq 4 ]]; then
    bowtie_input="-x $3 -1 $1 -2 $2"
    samtools_output="-o $4"   
else
	echo "Incorrect number of parameters"
	exit 2
fi
program_parameters="$program_parameters $bowtie_input"


#Run bowtie
eval "$program $program_parameters | samtools view -@ 20  -bSF4 --threads $SLURM_CPUS_ON_NODE - $samtools_output"


