#!/usr/bin/env bash

#SBATCH --mem=4G
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=4
#SBATCH --time=0:15:0


#Script for merging bam files using samtools samtools/1.12 on Compute Canada infrastructure using the local scratch directory
#If run as an sbatch job then the input files will be copied to the $SLURM_TMPDIR, merged, and the final result copied back
#If only one file is given to merge it will be copied directly and renamed
#Usage : sbatch .... localScratch_merge_bam_files.sh <merged_file> <bamfile1> [<bamfile2] ...


if [ "$#" -eq 2 ]; then
    mergedFile="$1"
    inputFile="$2"
    
    echo "Only one file to merge : copying $inputFile -> $mergedFile"
    cp $inputFile $mergedFile
else
    if [ -n "$SLURM_CPUS_ON_NODE" ]; then
        threads=$SLURM_CPUS_ON_NODE
    else
        threads=1
    fi

    if [ -n "$SLURM_TMPDIR" ]; then
        mergedFile="$1"
        scratchMergedFile="$SLURM_TMPDIR"/"$(basename $mergedFile)"
        scratchInputFiles=$(echo "${@:2}" | tr ' ' '\n' | xargs -I {} bash -c "echo $SLURM_TMPDIR/"'$(basename {})' | tr '\n' ' ' )
        
        echo "Copying input files ${@:2} ->  $SLURM_TMPDIR"
        echo "${@:2}" | tr ' ' '\n' | parallel --will-cite --eta -j "$threads" cp {} "$SLURM_TMPDIR"
        
        
        echo ./merge_bam_files.sh "$scratchMergedFile" "$scratchInputFiles"
        sleep 5
        
        echo "Copying $scratchMergedFile -> $mergedFile"
        cp "$scratchMergedFile" "$mergedFile"
                
        
    else
        mergedFile="$1"
        inputFiles="${@:2}"
        
        ./merge_bam_files.sh "$mergedFile" "$inputFiles"
    fi
fi

