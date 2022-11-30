#!/bin/bash

##################################################################
#
# Michael Mansfield, 2021
#
# This script pulls a number of singularity images which are
# used to execute all of the downstream commands needed to
# reproduce the analysis of the TeNT variants.
#
# Prerequisites:
#      -singularity installed and accessible on the PATH
# Usage:
#      ./setup.sh
#
#
##################################################################

case "$1" in
	-h|--h|-help|--help)
		echo
		echo "This script sets up a reproducible virtualized environment."
		echo "    Usage:     setup.sh [singularity]"
		echo "               If singularity is selected, please specify a singularity"
		echo "               cache directory. If none is specified, ./tools is used."
		exit 1
	;;
	*)
		METHOD="${}"
esac

# Find the directory of the setup.sh script, which lets you move
# between other directories easily.
SCRIPTDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
LIBDIR=$(echo ../lib)
TOOLDIR=$(echo ../tools)

cd "${LIBDIR}"
IMAGEFILE="images.txt"

printf "Pulling containers...\n\nTool\tVersion\tSource\n"
declare -A IMAGES
grep -v "#" "${IMAGEFILE}" | cut -d '"' -f 2 | while read LINE;
do
	#TOOL=$(echo "${LINE}" | cut -d "/" -f 3 | cut -d ":" -f 1)
	LINK=$(echo "${LINE}" | cut -d "\"" -f 2)
	TOOL=$(echo "${LINE}" | awk -F '/' '{print $NF}' | cut -d ':' -f 1)
	VERSION=$(echo "${LINE}" | awk -F '/' '{print $NF}' | cut -d ':' -f 2)
	printf "$TOOL\t$VERSION\t$LINK\n"
	IMAGES["${TOOL}"]+="${LINK}"

	# Pull relevant containers.
	for IMAGE in "${IMAGES[@]}";
	do
		singularity pull --dir "${TOOLDIR}" docker://"${IMAGE}"
	done
done
