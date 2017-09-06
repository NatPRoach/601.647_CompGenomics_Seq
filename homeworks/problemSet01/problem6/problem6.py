#!/usr/bin/env python

##########################################################################################################
#-----------------------------------------SCRIPT DESCRIPTION---------------------------------------------#
##########################################################################################################
# Script for global allignment of two arbitrary protein sequences specified in fasta format weighted by
# BLOSUM62 matrix
# Code written by Nathan P. Roach for 601.647 Computational Genomics: Sequences
# Contact info:
# Email 1: nroach2@jhu.edu
# Email 2: natproach@gmail.com
import sys
import string
import numpy as NUMPY
from operator import itemgetter
def scorePair62(residue1,residue2):
# Wrapper function for BLOSUM62 matrix scoring
    index = {'A':0, 'C':1, 'D':2, 'E':3, 'F':4, 'G':5, 'H':6, 'I':7, 'K':8, 'L':9, 'M':10, 'N':11, 'P':12,
        'Q':13, 'R':14, 'S':15, 'T':16, 'V':17, 'W':18, 'Y':19}
    blosum62Matrix=  []
#                           A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    blosum62Matrix.append([ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2])#A
    blosum62Matrix.append([ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2])#C
    blosum62Matrix.append([-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3])#D
    blosum62Matrix.append([-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2])#E
    blosum62Matrix.append([-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3])#F
    blosum62Matrix.append([ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3])#G
    blosum62Matrix.append([-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2])#H
    blosum62Matrix.append([-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1])#I
    blosum62Matrix.append([-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2])#K
    blosum62Matrix.append([-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1])#L
    blosum62Matrix.append([-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1])#M
    blosum62Matrix.append([-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2])#N
    blosum62Matrix.append([-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3])#P
    blosum62Matrix.append([-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1])#Q
    blosum62Matrix.append([-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2])#R
    blosum62Matrix.append([ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2])#S
    blosum62Matrix.append([ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2])#T
    blosum62Matrix.append([ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1])#V
    blosum62Matrix.append([-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2])#W
    blosum62Matrix.append([-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7])#Y
#                           A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    return blosum62Matrix[index[residue1]][index[residue2]]

def affineGapPenalty(a,b,L):
    return a + b*(L-1)

def dictatedPenaltyFunction(L):
    return affineGapPenalty(11,1,L)

def fastaReader(fasta_file):
    sequenceNames = []
    sequences = []
    sequenceCount = 0
    for line in fasta_file:
        print line
        if line[0] == '>':
            sequenceNames.append(string.rstrip(string.lstrip(line,'>'),'\n'))
            sequences.append('')
            sequenceCount += 1
        else:
            sequences[sequenceCount-1] = sequences[sequenceCount-1] + string.rstrip(line,'\n')
    return sequenceNames, sequences

def gotoh(a, b, scoreFunction, gapOpenPenalty, gapMaintainPenalty):
#generic implementation of the smith waterman allignment algorithm for strings a & b, using the score as dictated
#by score function.

    a_axisLen = len(a) + 1
    b_axisLen = len(b) + 1
##### Create matrix and initialize first row and collumn
    allignmentMatrix = NUMPY.zeros((a_axisLen,b_axisLen), dtype=int)
    aGapMatrix = NUMPY.zeros((a_axisLen,b_axisLen), dtype=int)
    bGapMatrix = NUMPY.zeros((a_axisLen,b_axisLen), dtype=int)
    directionMatrix = NUMPY.zeros((a_axisLen,b_axisLen), dtype=int) #direction matrix 0 = diagnal movement, 1 = a axis movement, 2 = b axis movement
    allignmentMatrix[0,0] = 0
    maxValue = NUMPY.iinfo(aGapMatrix.dtype).max
    aGapMatrix[0,0] = maxValue
    bGapMatrix[0,0] = maxValue
    for i in range(1,a_axisLen):
        allignmentMatrix[i,0] = gapOpenPenalty + gapMaintainPenalty*(i-1)
        bGapMatrix[i,0] = allignmentMatrix[i,0]
    for j in range(1,b_axisLen):
        allignmentMatrix[0,j] = gapOpenPenalty + gapMaintainPenalty*(j-1)
        aGapMatrix[0,j] = allignmentMatrix[0,j]
##### Push out matrix
    globalMax = 0
    globalMaxCoordinates = (0,0)
    for i in range(1,a_axisLen):
        for j in range(1,b_axisLen):
            aGapMatrix[i,j] = max(allignmentMatrix[i-1,j] + gapOpenPenalty, aGapMatrix[i-1,j] + gapMaintainPenalty)
            bGapMatrix[i,j] = max(allignmentMatrix[i,j-1] + gapOpenPenalty, bGapMatrix[i,j-1] + gapMaintainPenalty)
            (directionMatrix[i,j], allignmentMatrix[i,j]) = max(enumerate((allignmentMatrix[i-1,j-1]+scoreFunction(a[i-1],b[j-1]),
                aGapMatrix[i,j], bGapMatrix[i,j])),key=itemgetter(1))
            if allignmentMatrix[i,j] > globalMax:
                globalMax = allignmentMatrix[i,j]
                globalMaxCoordinates = (i,j)
##### Walk back
    a_prime = []
    b_prime = []
    i = a_axisLen-1
    j = b_axisLen-1
    
    print allignmentMatrix[i,j]

    while i > 0 and j > 0:
        direction = directionMatrix[i,j]
        if direction == 0:
            a_prime.append(a[i-1])
            b_prime.append(b[j-1])
            i -= 1
            j -= 1
        elif direction == 1:
            a_prime.append(a[i-1])
            b_prime.append('-')
            i -= 1
        elif direction == 2:
            a_prime.append('-')
            b_prime.append(b[j-1])
            j -= 1
    
    a_prime.reverse()
    b_prime.reverse()
    print ''.join(a_prime)
    print ''.join(b_prime)

input_file = open(sys.argv[1])
sequenceNames, sequences = fastaReader(input_file)
gotoh(sequences[0],sequences[1],scorePair62,-11,-1)