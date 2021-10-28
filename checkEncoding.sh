#!/usr/bin/env bash

#Small script for predicting the encoding of fastq files
#Modified from https://bioinformaticsworkbook.org/introduction/fastqquality-score-encoding.htm
#Usage : checkEncoding.sh <fastq_file>

head -n 10000 "$1" |\
  awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 -v | \
  awk 'BEGIN{min=100;max=0;} \
      {for(i=1;i<=NF;i++) \
          {if($i>max) max=$i; \
               if($i<min) min=$i;}}END \
          {if(max<=74 && min<59) \
                     print "Phred+33"; \
           else \
           if(max>73 && min>=64) \
                     print "Phred+64"; \
           else \
           if(min>=59 && min<64 && max>73) \
                     print "Solexa+64"; else print "Unknown score encoding!", min," ", max;}'
