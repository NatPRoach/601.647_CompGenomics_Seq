#!/usr/bin/env python

##########################################################################################################
#-----------------------------------------SCRIPT DESCRIPTION---------------------------------------------#
##########################################################################################################
# Script for transcribing an arbitrary DNA sequence.
# Code written by Nathan P. Roach for 601.647 Computational Genomics: Sequences
# Contact info:
# Email 1: nroach2@jhu.edu
# Email 2: natproach@gmail.com
import sys
import string
input_file = open(sys.argv[1])
input_file_string = input_file.read()
print string.rstrip(string.replace(input_file_string,'T','U'),'\n')