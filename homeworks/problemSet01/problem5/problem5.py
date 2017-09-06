#!/usr/bin/env python

##########################################################################################################
#-----------------------------------------SCRIPT DESCRIPTION---------------------------------------------#
##########################################################################################################
# Script for translating arbitrary RNA sequence.
# Code written by Nathan P. Roach for 601.647 Computational Genomics: Sequences
# Contact info:
# Email 1: nroach2@jhu.edu
# Email 2: natproach@gmail.com
import sys
import string
def proteinTranslate(mRNA_string):
# Function wrapper for translating mRNA. 
    codonNum = len(mRNA_string)/3
    proteinString = ['']*codonNum
    for i in range(codonNum):
        proteinString[i] = codonCall(mRNA_string[i*3:i*3+3])
    return ''.join(proteinString)

def codonCall(codon):
# Function wrapper for calling amino acid for a given codon
    codonTable = {
    'UUU':'F','CUU':'L','AUU':'I','GUU':'V',
    'UUC':'F','CUC':'L','AUC':'I','GUC':'V',
    'UUA':'L','CUA':'L','AUA':'I','GUA':'V',
    'UUG':'L','CUG':'L','AUG':'M','GUG':'V',
    'UCU':'S','CCU':'P','ACU':'T','GCU':'A',
    'UCC':'S','CCC':'P','ACC':'T','GCC':'A',
    'UCA':'S','CCA':'P','ACA':'T','GCA':'A',
    'UCG':'S','CCG':'P','ACG':'T','GCG':'A',
    'UAU':'Y','CAU':'H','AAU':'N','GAU':'D',
    'UAC':'Y','CAC':'H','AAC':'N','GAC':'D',
    'UAA':'.','CAA':'Q','AAA':'K','GAA':'E',
    'UAG':'.','CAG':'Q','AAG':'K','GAG':'E',
    'UGU':'C','CGU':'R','AGU':'S','GGU':'G',
    'UGC':'C','CGC':'R','AGC':'S','GGC':'G',
    'UGA':'.','CGA':'R','AGA':'R','GGA':'G',
    'UGG':'W','CGG':'R','AGG':'R','GGG':'G'}
    
    return codonTable[codon]

input_file = open(sys.argv[1])
input_file_string = input_file.read()

print string.rstrip(proteinTranslate(input_file_string),'.')