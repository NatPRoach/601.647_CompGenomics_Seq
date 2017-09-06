#!/usr/bin/env python

##########################################################################################################
#-----------------------------------------SCRIPT DESCRIPTION---------------------------------------------#
##########################################################################################################
# Script for counting number of each nucleotide within an arbitrary DNA sequence.
# Code written by Nathan P. Roach for 601.647 Computational Genomics: Sequences
# Contact info:
# Email 1: nroach2@jhu.edu
# Email 2: natproach@gmail.com
import sys
input_file = open(sys.argv[1])
input_file_string = input_file.read()
nucleotideCount = [0, 0, 0, 0] # A C G T in that order
for char in input_file_string:
    if char == 'A':
        nucleotideCount[0] += 1
    elif char == 'C':
        nucleotideCount[1] += 1
    elif char == 'G':
        nucleotideCount[2] += 1
    elif char == 'T':
        nucleotideCount[3] += 1

print "%d %d %d %d" %(nucleotideCount[0],nucleotideCount[1],nucleotideCount[2],nucleotideCount[3])