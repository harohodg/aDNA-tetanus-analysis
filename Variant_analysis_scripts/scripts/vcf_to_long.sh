find merged singles -iname '*.vcf' | while read i; do b=$(basename $i | awk -F '.octopus' '{print $1}'); grep -v '^#' $i | while read v; do printf "$b\t$v\n"; done; done
