#!/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import sys
import os

def check_positive_integer(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is not 0 or a positive integer" % value)
    return ivalue

parser = argparse.ArgumentParser(description='Don\'t use this script for anything. It\'s really just for one very specific purpose, which is to mask low/no-coverage positions from a reference FASTA file before calling consensus sequences from a VCF file.', formatter_class=RawTextHelpFormatter)
parser.add_argument('-f','--fastaFile', help='Reference FASTA file to be masked. It assumes the FASTA file has one entry, for now.\nDefault:none', required=True)
parser.add_argument('-d','--depthFile', help='Output of samtools depth, with reads aligned to the reference FASTA file.\nDefault: none', required=True)
parser.add_argument('-o','--outFile', help='Output file with masked positions.\nDefault: [input fasta].masked.fa', required=False)
parser.add_argument('--lowCoverageCutoff', help='Cutoff for low-coverage masking. Positions with coverage equal or less than this value are masked with the character determined by lowCoverageChar.\nDefault: 5', required=False, default=5, type=check_positive_integer)
parser.add_argument('--lowCoverageChar', help='Character to mask low-coverage positions with.\nDefault: N. Special option: \'lower\', which makes low-coverage positions lowercase.', required=False, default='N')
parser.add_argument('--noCoverageCutoff', help='Positions with coverage equal or less than this value are masked with the character determined by noCoverageChar.\nDefault: 0', required=False, default=0, type=check_positive_integer)
parser.add_argument('--noCoverageChar', help='Character to mask low-coverage positions with.\nDefault: -.', required=False, default='-')
parser.add_argument('--region', help='Region of the reference FASTA to mask. Useful when the depth file represents a subset of the reference FASTA.\nDefault: none. Format: \'start-end\'', default=0, required=False)
parser.add_argument('--maskRefLowercase', help='Mask reference positions that pass the filter as lowercase.\nDefault: False', required=False, action='store_true')
args = parser.parse_args()

if args.outFile:
    out_fasta = args.outFile
else:
    out_fasta = os.path.basename(args.fastaFile) + '.masked.fa'

if args.region:
    try:
        region = args.region.replace('"', '').replace("'",'')
        region_start = int(region.split('-')[0])
        region_end = int(region.split('-')[1])
    except:
        print('Error. Are you sure that you specified the region correctly? The format should be: \'start-end\' or similar.')
        sys.exit(1)

if args.lowCoverageCutoff <= args.noCoverageCutoff:
    print('Error.-- noCoverageCutoff must be less than --lowCoverageCutoff.')
    sys.exit(1)

if len(args.lowCoverageChar) > 1:
    if args.lowCoverageChar == 'lower':
        mask_low_cov_lower = True
    else:
        print('Error. --lowCoverageChar must be a single character or the special argument \'lower\'.')
        sys.exit(1)

# Dictionary structure: {fasta_header:seq}
fasta_dict = {}
with open(args.fastaFile) as fasta_handler:
    for line in fasta_handler:
        line = line.strip()
        if line.startswith('>'):
            header = line.split('>')[1]
            # This line splits the header.
            # Note that these splits will probably not always work.
            header = header.split(' ')[0].split(':')[0]
            fasta_dict[header] = ''
        else:
            fasta_dict[header] += line

# Dictionary structure: {sequence_name:[position, depth]}
depth_dict = {}
with open(args.depthFile) as depth_handler:
    counter = 0
    for line in depth_handler:
        line = line.strip()
        positionals = line.split('\t')
        if positionals[0] not in depth_dict:
            depth_dict[positionals[0]] = [[positionals[1], positionals[2]]]
        else:
            depth_dict[positionals[0]].append([positionals[1], positionals[2]])

with open(out_fasta, 'w') as out_handler:
    for fasta_header in fasta_dict:
        sequence = fasta_dict[fasta_header]
        try:
            depth_info = depth_dict[fasta_header]
        except:
            print('Error. It looks like the fasta header %s has no equivalent in the depth file?' % fasta_header)
            sys.exit(1)

        # Construct a sequence to return, based on what the coverage value for that position is.
        # Does the subsetting of both the FASTA file and the depth
        # file here.
        return_sequence = ''
        depth_values = [int(x[1]) for x in depth_info]
        if region:
            #subsequence = sequence[region_start-1 : region_end]
            depth_values = depth_values[region_start-1 : region_end-1]
            position = 0
            for depth in depth_values:
                if depth <= args.lowCoverageCutoff:
                    if depth <= args.noCoverageCutoff:
                        letter = args.noCoverageChar
                    else:
                        try:
                            if mask_low_cov_lower:
                                letter = sequence[position].lower()
                        except:
                            letter = args.lowCoverageChar
                else:
                    if args.maskRefLowercase:
                        letter = sequence[position].lower()
                    else:
                        letter = sequence[position]
                return_sequence += letter
                position += 1
        # If no subsetting is necessary, ensure that the depth dictionary and FASTA sequence are the same length before anything.
        else:
            if len(depth_values) != len(sequence):
                print('Error. It looks like the number of depths and length of the FASTA sequence %s are not equal.' % fasta_header)
                sys.exit(1)
            else:
                position = 0
                for depth in depth_values:
                    if depth <= args.lowCoverageCutoff:
                        if depth <= args.noCoverageCutoff:
                            letter = args.noCoverageChar
                        else:
                            try:
                                if mask_low_cov_lower:
                                    letter = sequence[position].lower()
                            except:
                                letter = args.lowCoverageChar
                    else:
                        if args.maskRefLowercase:
                            letter = sequence[position].lower()
                        else:
                            letter = sequence[position]
                    return_sequence += letter
                    position += 1
        # Format the output sequence to have nice newlines every X letters.
        #return_sequence = ''.join(return_sequence[i:i+64] + '\n' for i in range(0, len(return_sequence), 60))
        out_handler.write('>%s\n%s\n' % (fasta_header, return_sequence))
