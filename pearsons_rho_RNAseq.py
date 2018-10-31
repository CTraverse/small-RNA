#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 15:20:37 2018

@author: chuck
"""

"""
This script will calculate Person's R for each RNAseq library to determien how close each replicate is to one another.
It appears that there's a systematic bias in replicate 1, so we will test this by looking at genome-wide expression correlation
"""


from scipy.stats.stats import pearsonr   
from sklearn import preprocessing
import numpy as np
import pandas as pd

#import the reads and the metadata for each RNAseq library
read_counts = pd.read_csv('/home/chuck/Documents/RNAseq/coding_norm/coding_gene_sums_norm.txt', sep='\t')
read_information = pd.read_csv('/home/chuck/Documents/RNAseq/RNAseq_sample_info/ltee_clone_names.txt', sep='\t')

#Keep only the read counts
read_counts = read_counts.drop(['features'], 1)

#Get a list of the bacterial clones corresponding to each RNAseq library
clones = list(set(read_information["clone"]))

#Create an empty dictionary of dictionaries to call the correct replicate for each clone by its clone name and replicate number
clone_dict = {}
for clone in clones:
    clone_dict[clone] = {}
    for i in range(3):
        clone_dict[clone][i+1] = ""
        
#Populate the dict of dicts with the unique library identifier
for library in read_information.iterrows():
    clone_dict[library[1]["clone"]][int(library[1]["replicate"])] = library[1]["sample"]
    
#Standardize the reads counts so each library is comparable to one another
standard = pd.DataFrame(preprocessing.scale(read_counts), columns = read_counts.columns)

r_write = open('/home/chuck/Documents/RNAseq/coding_norm/pearsonr_library_replicates.txt', 'w')

#Calculate Pearson's R to find how well correlated the library replicates are to one another
for clone in clones:
    first, second, third = read_counts[clone_dict[clone][1]], read_counts[clone_dict[clone][2]], read_counts[clone_dict[clone][3]]
    #Calculate for each combination of the two replicates
    oneVStwo = pearsonr(first,second)
    oneVSthree = pearsonr(first,third)
    twoVSthree = pearsonr(second,third)
    
    #Write the Pearson's Rs to a file
    r_write.write("%s\t%s\t%s\t%s" % (clone, oneVStwo[0], oneVSthree[0], twoVSthree[0]))





