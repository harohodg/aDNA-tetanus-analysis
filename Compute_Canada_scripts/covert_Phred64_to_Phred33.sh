#!/usr/bin/env bash

#Usage <input file> <output file>
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

$DIR/seqtk_job.sh 'seq -VQ64' "$1" "$2"
