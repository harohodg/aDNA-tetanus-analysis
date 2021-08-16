#!/usr/bin/env bash

#Script to calculate trimmomatic job time in minutes given input file size
#$1 = min time in minutes
#$2 = max time in minutes
#$3 = max file(s) size
#$4 = current file(s) size

echo "scale=3;a=($1+($4/$3)*($2 - $1));scale=0;a/1+1" | bc
