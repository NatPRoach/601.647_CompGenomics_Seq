#!/usr/bin/env python

##########################################################################################################
#-----------------------------------------SCRIPT DESCRIPTION---------------------------------------------#
##########################################################################################################
# Script for creating the reverse compliment of arbitrary DNA sequence.
# Code written by Nathan P. Roach for 601.647 Computational Genomics: Sequences
# Contact info:
# Email 1: nroach2@jhu.edu
# Email 2: natproach@gmail.com
import sys
import string
input_file = open(sys.argv[1])
input_file_string = string.rstrip(input_file.read(), ' \n')
seq_length = len(input_file_string)
rcomp_list = ['']*seq_length
for i, char in enumerate(input_file_string,start=1):
    if char == 'A':
        rcomp_list[seq_length-i] = 'T'
    elif char == 'C':
        rcomp_list[seq_length-i] = 'G'
    elif char == 'G':
        rcomp_list[seq_length-i] = 'C'
    elif char == 'T':
        rcomp_list[seq_length-i] = 'A'

print ''.join(rcomp_list)