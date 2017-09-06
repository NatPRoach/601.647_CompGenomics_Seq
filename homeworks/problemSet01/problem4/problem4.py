#!/usr/bin/env python

##########################################################################################################
#-----------------------------------------SCRIPT DESCRIPTION---------------------------------------------#
##########################################################################################################
# Script for computing the hamming distance of two arbitrary DNA sequences of equal length
# Code written by Nathan P. Roach for 601.647 Computational Genomics: Sequences
# Contact info:
# Email 1: nroach2@jhu.edu
# Email 2: natproach@gmail.com
import sys
input_file = open(sys.argv[1])
input_file_sequence1 = input_file.readline()
input_file_sequence2 = input_file.readline()

hammingDistance = 0
for i in range(len(input_file_sequence1)):
    if input_file_sequence1[i] != input_file_sequence2[i]:
        hammingDistance += 1

print hammingDistance