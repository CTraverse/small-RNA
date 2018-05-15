#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 15:20:37 2018

@author: chuck
"""

from scipy.stats.stats import pearsonr   
from sklearn import preprocessing
import numpy as np
import pandas as pd

read_counts = pd.read_csv('/home/chuck/Documents/RNAseq/coding_norm/coding_gene_sums_norm.txt', sep='\t')
read_information = pd.read_csv('/home/chuck/Documents/RNAseq/RNAseq_sample_info/ltee_clone_names.txt', sep='\t')

read_counts = read_counts.drop(['features'], 1)

clones = list(set(read_information["clone"]))

clone_dict = {}
for clone in clones:
    clone_dict[clone] = {}
    for i in range(3):
        clone_dict[clone][i+1] = ""

for library in read_information.iterrows():
    clone_dict[library[1]["clone"]][int(library[1]["replicate"])] = library[1]["sample"]
    

#print(read_counts)
#print(clone_dict)
standard = pd.DataFrame(preprocessing.scale(read_counts), columns = read_counts.columns)

for clone in clones:
    first, second, third = read_counts[clone_dict[clone][1]], read_counts[clone_dict[clone][2]], read_counts[clone_dict[clone][3]]
    oneVStwo = pearsonr(first,second)
    oneVSthree = pearsonr(first,third)
    twoVSthree = pearsonr(second,third)

    print("%s\t%s\t%s\t%s" % (clone, oneVStwo[0], oneVSthree[0], twoVSthree[0]))





