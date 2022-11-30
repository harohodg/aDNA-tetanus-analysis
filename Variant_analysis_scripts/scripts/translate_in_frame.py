#!/bin/python


import argparse
from argparse import RawTextHelpFormatter
import sys
import os

parser = argparse.ArgumentParser(description='Don\'t use this script for anything! The purpose of this script is to make a perfect in-frame translation of a nucleotide MSA.', formatter_class=RawTextHelpFormatter)
parser.add_argument('-f','--fastaFile', help='Multiple alignment of reference nucleotide sequences and consensus variants.', required=True)
parser.add_argument('-of','--outputFasta', help='Output FASTA file for in-frame translations.\nDefault: in_frame_translation.fasta', required=False, default='translated_in_frame.fasta')
parser.add_argument('-ot','--outputTSV', help='Output file with non-synonymous and synonymous subtitutions marked.\nDefault:none', required=False, default='translated_in_frame.substitutions.tsv')
parser.add_argument('--reportStopsAsNonsynonymous', help='Whether or not to report stop codons as being non-synonymous substitutions.\nDefault: False', required=False, action='store_true')
args = parser.parse_args()

# Dictionary structure: {fasta_header:seq}
fasta_dict = {}
with open(args.fastaFile) as fasta_handler:
    for line in fasta_handler:
        line = line.strip()
        if line.startswith('>'):
            header = line.split('>')[1]
            fasta_dict[header] = ''
        else:
            fasta_dict[header] += line

standard_genetic_code = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

def split_into_codons(seq):
    subseq = []
    for i in range(0, len(seq), 3):
        subseq.append(seq[i:i+3].upper())
    return(subseq)

def translate_nt_with_missing(triplet):
    if triplet in standard_genetic_code:
        aa = standard_genetic_code[triplet]
    elif triplet == "---":
        aa = "-"
    else:
        aa = "X"
    return(aa)

def find_substitutions(reference_sequence, query_sequence):
    reference_codons = split_into_codons(reference_sequence)
    query_codons = split_into_codons(query_sequence)

    substitution_dictionary = {}
    for codon in range(0, len(reference_codons)):
        reference_codon = reference_codons[codon]
        reference_translation = translate_nt_with_missing(reference_codon)
        query_codon = query_codons[codon]
        query_translation = translate_nt_with_missing(query_codon)
        if args.reportStopsAsNonsynonymous:
            if query_codon != "---" and query_translation != "X":
                if reference_codon != query_codon:
                    if reference_translation != query_translation:
                        substitution_dictionary[codon] = {
                            "Codon" : codon,
                            "Ref codon" : reference_codons[codon],
                            "Ref aa" : reference_translation,
                            "Query codon" : query_codons[codon],
                            "Query aa" : query_translation,
                            "Subs. type" : "Non-synonymous"
                            }
                    else:
                        substitution_dictionary[codon] = {
                            "Codon" : codon,
                            "Ref codon" : reference_codons[codon],
                            "Ref aa" : reference_translation,
                            "Query codon" : query_codons[codon],
                            "Query aa" : query_translation,
                            "Subs. type" : "Synonymous"
                            }
        else:
            if query_codon != "---" and query_translation != "X" and query_translation != "*":
                if reference_codon != query_codon:
                    if reference_translation != query_translation:
                        substitution_dictionary[codon] = {
                            "Codon" : codon,
                            "Ref codon" : reference_codons[codon],
                            "Ref aa" : reference_translation,
                            "Query codon" : query_codons[codon],
                            "Query aa" : query_translation,
                            "Subs. type" : "Non-synonymous"
                            }
                    else:
                        substitution_dictionary[codon] = {
                            "Codon" : codon,
                            "Ref codon" : reference_codons[codon],
                            "Ref aa" : reference_translation,
                            "Query codon" : query_codons[codon],
                            "Query aa" : query_translation,
                            "Subs. type" : "Synonymous"
                            }
    return(substitution_dictionary)

translated_dictionary = {}
for fasta_header in fasta_dict:
    split_seq = split_into_codons(fasta_dict[fasta_header])
    translated_seq = ""
    for triplet in split_seq:
        translated_seq += translate_nt_with_missing(triplet)
    translated_dictionary[fasta_header + "_translated"] = translated_seq

with open(args.outputFasta, "w") as o:
    for seq in translated_dictionary:
        o.write(">%s\n%s\n" % (seq, translated_dictionary[seq]))

with open(args.outputTSV, "w") as o:
    o.write("Sample\tCodon_number\tE88_codon\tE88_AA\tSample_codon\tSample_AA\tSubstitution_type\n")
    for seq in fasta_dict:
        if seq != "E88_Reference":
            substitutions = find_substitutions(fasta_dict["E88_Reference"], fasta_dict[seq])
            for substitution in substitutions:
                o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    seq,
                    substitutions[substitution]['Codon'],
                    substitutions[substitution]['Ref codon'],
                    substitutions[substitution]['Ref aa'],
                    substitutions[substitution]['Query codon'],
                    substitutions[substitution]['Query aa'],
                    substitutions[substitution]['Subs. type']
                    )
                )
