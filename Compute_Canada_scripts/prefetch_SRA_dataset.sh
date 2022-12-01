#!/usr/bin/env bash

#Script for prefetching a single SRA dataset
#$1 the SRA ID
#$2 the target folder

module load nixpkgs/16.09 sra-toolkit/2.9.6

dataSize=$( vdb-dump --info "$1" | grep "size" | cut -f2 -d':' | tr -d ', ')
maxSize=$(expr "$dataSize" \* 2)
prefetch --max-size "$maxSize" -O "$2" "$1"
