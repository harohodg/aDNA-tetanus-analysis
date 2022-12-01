#!/usr/bin/env bash


#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --job-name=parallel_mpileup_jobs
#SBATCH --output=%x-%j.out



#Script for running samtools/1.15.1 mpileup on Compute Canada infrastructure
#If -d flag is used the program will print what what would have been run to stdout
#If -c compression_level flag is used the output will be gz compressed with the same level
#If compression level is set to zero no compression will be applied
#If -t num_threads is given, this many threads will be used with pigz
#Parameters are assumed to be output file, bam file, optional reference file
# !!! DOES NOT CHECK FOR EXISTING OUTPUT FILE !!!
#Author : Harold Hodgins <hhodgins@uwaterloo.ca>

#Usage : [sbatch [--account=?] [--time=?] [--job-name=?]] mpileup_job.sh [-d]  [-t threads (default 1)] [-c compression_level 0-9 (0=None, default 9)] <bamFile> <outputFile>  [any parameters to pass to samtools mpileup]


#Version History 
#Version 1.0 : September 6, 2022
#   Functional script with minimal error checking


VERSION='1.0.0'
DEFAULT_COMPRESSION_THREADS='1' #How many threads do we use for compression
DEFAULT_COMPRESSION_LEVEL=9

>&2 echo "mpileup_job.sh version $VERSION"

# Echo usage if something isn't right.
usage() { 
    echo "Usage: $0 [-d] [-t threads (default ${DEFAULT_COMPRESSION_THREADS})] [-c compression_level 0-9 (0=None, default ${DEFAULT_COMPRESSION_LEVEL})] <bamFile> <outputFile>  [any parameters to pass to samtools mpileup] " 1>&2; exit 1; 
}


while getopts ":dt:c:" o; do
    case "${o}" in
        d)  
            debug=1
            ;;
        c)
            compressionLevel="$OPTARG"
            re='^[0-9]$'
            if ! [[ $compressionLevel =~ $re ]] ; then
               echo "error: compression level is not between 0-9" >&2; exit 1
            fi
            ;;
        t)
            compressionThreads="$OPTARG"
            re='^[0-9]+$'
            if ! [[ $compressionThreads =~ $re ]] ; then
               echo "error: compression threads is not a number" >&2; exit 1
            fi
            ;;
        \?)
            echo "ERROR: Invalid option -$OPTARG"
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done
shift $((OPTIND-1))



if [[ "$#" -lt 2 ]]; then
    echo 'Incorrect number of arguments.'
    usage
else
    bamFile="$1"
    outputFile="$2"
    mpileupParameters="${@:3}"
fi


compressionLevel="${compressionLevel:-$DEFAULT_COMPRESSION_LEVEL}"
compressionThreads="${compressionThreads:-$DEFAULT_COMPRESSION_THREADS}"
tmpOutputFile="${outputFile}.partial"


pigzParamaters="-${compressionLevel} --to-stdout --processes ${compressionThreads}"
mpileupParameters="${mpileupParameters} ${bamFile}"

mpileupCommand="samtools mpileup ${mpileupParameters}"
addQuotations=$'sed "s|\\t|\'\\t\'|g" | awk -v quote="\'" \'{print quote $0 quote}\''
if [[ $compressionLevel -ne 0 ]];then
    mpileupCommand="${mpileupCommand} |  ${addQuotations} | pigz ${pigzParamaters} > ${tmpOutputFile}"
else
    mpileupCommand="${mpileupCommand} |  ${addQuotations}  > ${tmpOutputFile}"
fi
mpileupCommand="$mpileupCommand && mv $tmpOutputFile $outputFile"

jobSetupCommand="set -e;set -o pipefail; module load StdEnv/2020 samtools/1.15.1 && mkdir -p $(dirname $tmpOutputFile)"
jobCommand="$jobSetupCommand && $mpileupCommand"


if [[ -n $debug ]];then
    echo "$jobCommand"
else
    echo "Running samtools/1.15.1 mpileup with $jobCommand"
    eval "$jobCommand" && echo -n 'DONE'  || echo -n 'FAILED'
    echo " running samtools mpileup"
fi
