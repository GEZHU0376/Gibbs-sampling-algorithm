#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
this program contains list of modules can perform calculation functions. 
all module have been tested and successfully ran on vim and 
Spyder python 3.7 at local directory.

Created on Sat Apr 13 13:26:50 2019
@author: gezhu
"""
#imports 
import random
import numpy
import numpy as np
from scipy.stats import entropy
from math import log, e
import pandas as pd
from itertools import permutations
from itertools import combinations

# a. count
#Function to count the occurrence of each letter at every positions.
def Count(Motifs,p):
    #list of their counts at each position
    Count = {}

    # length of single motif
    k = len(Motifs[0]) 
    
    #get all symbols in motif:
    all_symbols = ""
    for motif in Motifs:
        for i in motif:
            if i not in all_symbols:
                all_symbols=all_symbols+i
    
    # Sort symbols alphabetically.
    all_symbols = list(all_symbols)
    all_symbols.sort()         
    all_symbols = ''.join(all_symbols)
        
    
    #check value if is pseudo-count
    if p==0:
        for letter in all_symbols:
            Count[letter]=[0.0] * k
    else:
        for letter in all_symbols:
            Count[letter]=[float(p)] * k
            
    #Append dictionary for each symbol encountered in the process
    for motif in Motifs:
        for j in range(len(motif)):
            symbol = motif[j]
            Count[symbol][j] += 1
            
    return Count

#print(Count(["ACGT","AAAC","AGTT","CCGT"],1))

#b. profile
def profile(string_list,p=0):
    t = len(string_list)
    k = len(string_list[0]) 
    profile = Count(string_list,p)
    
    for symbol in profile:
        for i in range(k):
            profile[symbol][i] = profile[symbol][i]/t
            
    return profile

string_list = ["AABDCA","CCCACB","CACCDB","AACCBB","AACCDD"]   
#print (profile(string_list,2))

#c. substring
#Function to generate all the possible substrings from input
"""Function to return the possible substrings of length k"""
def substrings(string,k):
    substringList=list()
    for i in range(k,len(string)+1):
        substringList.append(string[i-k:i])
    return substringList
#    output=list(permutations(string,k)) 
#    #produces a list containing tuples in specific order
#    outputlist=list()
#    #list that store final output
#    for x in output:
#        #Convert tuples in list to string
#        stringoutput=  ''.join(x)
#        outputlist.append(stringoutput)
#    return outputlist
    
#example
#print (substrings("ACGT",2))

#d. prob_motif
#Function to calculate the probability of input string/motif
def prob_motif (motif, profile):
    prob = 1.0
    
    for i in range(len(motif)):
        symbol = motif[i]
        if symbol not in profile:
            prob *= 0
        else:
            prob *= profile[symbol][i]

    return prob
    
#print( prob_motif("abc", {'a':0.25, 'b':0.25, 'c': 0.25, 'd': 0.25}) )

#e. select_new
#Function to select a random value from motifs based on weighted probability
#Function to normalize probailities if not equal to one
def normalize(probabilities):
    factor = 1 / sum(probabilities)
    return [factor * p for p in probabilities]

#Function to find the interval
def interval(x, partition):
    for i in range(0, len(partition)):
        if x < partition[i]: 
            return i-1
    return -1

def select_new(motifs,probs):
    #check if sum of probability is equal to one. If not normalize
    if sum(probs) != 1:
        probs=normalize(probs)
    
    x = numpy.random.random()
    cum_weights = [0] + list(numpy.cumsum(probs))
    index = interval(x, cum_weights)
    return motifs[index]

#Example
#print(select_new(["AACGTA","CCCGTT","CACCTT","TTCCGG","TTCCGG"],[1/12, 1/6, 3/6, 1/6, 1/6]))
#output -> CACCTT

#f. entropy:
def entropy(labels, base=None):
  value,counts = np.unique(labels, return_counts=True)
  norm_counts = counts / counts.sum()
  base = e if base is None else base
  return -(norm_counts * np.log(norm_counts)/np.log(base)).sum()

z = entropy(['ACGGT', 'TCGGA',  'GGTT'])
#print(z)

# g. consensus count most frequent letter at each position of the list
def consensus(motifs):
    profileDict = {}
    motif_length = len(motifs[0])            
    for motif in motifs:       
        for i in range(len(motif)):     
            if motif[i] in profileDict:     
                profileDict[motif[i]][i]+=1
            else:                           
                profileDict[motif[i]] = [0]*motif_length
                profileDict[motif[i]][i]+=1
                
    consensus_str=""
    for i in range(motif_length):
        val = 0.0
        symbol = next(iter(profileDict.values()))
        for item in profileDict:
            if(profileDict[item][i] > val):
                val = profileDict[item][i]
                symbol = item
        consensus_str += symbol
            
    return consensus_str


#h. print_matrix(profile) 
# input profile of probability and output whole matrix 
def print_matrix(profile):
    for item in profile:
        print('%s ' %item, end="")
        for i in profile[item]:
            print('%.2f ' %i, end = "")
        print()
        
mymotifs = ['ACGTT','CCGTT', 'AGGTT', 'ACTTT', 'ACGAT','TTTTT']
my_profile = profile(mymotifs)
#print(my_profile)
#print_matrix(my_profile)
#print(consensus(mymotifs))

