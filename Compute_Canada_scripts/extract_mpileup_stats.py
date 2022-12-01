#!/usr/bin/env python3

#Small script for extracting various stats from mpileup results
#From a few test cases 10,000 rows per chunk uses < 512MB of RAM 


#Version History 
#Version 1.0 : September 23, 2022
#   Functional script with minimal error checking
#
#Version 2.0 : October 25, 2022
#   Removed phred filtering
#   Added mapped reads lengths
#
#Version 2.1 : November 19, 2022
#   Added mapped reads
PILEUP_COLUMNS = ('reference',
                  'position',
                  'reference_character',
                  'num_reads',
                  'mapped_characters',
                  'base_quality',
                  'mapping_quality',
                  'read_char_position',
                  'reads')
DATA_TYPES     = {'reference':'str', 
                  'position':'int',
                  'reference_character':'str',
                  'num_reads':'int',
                  'mapped_characters':'str',
                  'base_quality':'str',
                  'mapping_quality':'str',
                  'read_char_position':'str',
                  'reads':'str'}
COLUMNS_TO_INCLUDE = ('reference','position','num_reads','mapped_characters','base_quality', 'reads')        
DEFAULT_ROWS_PER_CHUNK = 1000
DEFAULT_PHRED_THRESHOLD=0

import sys
import json
import argparse
import pandas as pd
from collections import Counter


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--rows', type=int,    help='How many rows per chunk?', default=DEFAULT_ROWS_PER_CHUNK)
parser.add_argument('--output_file','-o',  help='Where do we write the json output to? (default = stdout)', type=argparse.FileType('w'), default=sys.stdout )
#parser.add_argument('--phred_score', '-p', type=int, help='What phred_score cuttoff are we using?', default=DEFAULT_PHRED_THRESHOLD)
parser.add_argument("input_file")

args  = parser.parse_args()
#assert args.phred_score >= 0


file_stats = {}
pileup_data = pd.read_csv( args.input_file, 
                           dtype=DATA_TYPES, 
                           header=None, 
                           compression='infer', 
                           sep='\t', 
                           names=PILEUP_COLUMNS, 
                           quotechar="'", 
                           keep_default_na=False, 
                           iterator=True, 
                           chunksize=args.rows, 
                           usecols=COLUMNS_TO_INCLUDE)
for chunk in pileup_data:
    for reference in chunk['reference'].unique():
        file_stats[reference]                 = file_stats.get(reference, {} )
        file_stats[reference]['phred_scores'] = file_stats[reference].get('phred_scores', Counter() )
        file_stats[reference]['coverage']     = file_stats[reference].get('coverage', Counter() )
        file_stats[reference]['read_lengths'] = file_stats[reference].get('read_lengths', Counter() )
        
        tmp_data = chunk.query('reference == "{}"'.format(reference)).copy()
        tmp_data['mapped_characters'] = tmp_data['mapped_characters'].apply(lambda s : s.replace(',','.').upper())
        tmp_data.apply(lambda row : file_stats[reference]['coverage'].update( {(0,0):1} ) if row.num_reads == 0 else
                ( file_stats[reference]['coverage'].update( {(str(row['mapped_characters']).count('.'),row.num_reads) :1}) ,
                file_stats[reference]['phred_scores'].update([ord(b)-33 for b in row.base_quality]),
                file_stats[reference]['read_lengths'].update(str(row.reads).split(','))
                ), axis=1 )

for reference in file_stats.keys():
    file_stats[reference]['reads'] = tuple( file_stats[reference]['read_lengths'].keys() )
    file_stats[reference]['read_lengths'] = Counter( file_stats[reference]['read_lengths'].values() )
    for data_label in ('phred_scores','coverage', 'read_lengths'):
        file_stats[reference][data_label] = tuple( file_stats[reference][data_label].items() )

json.dump( file_stats, args.output_file)
args.output_file.close()
