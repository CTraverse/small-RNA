#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 14:43:36 2018

@author: chuck
"""

"""
This script will take each RNAseq read counts file and generate three files: 
Sense reads, antisense reads, and a file that is a combination of both.
"""

import pandas as pd
import os

#Import list of files for each library
files = os.listdir("/home/chuck/Documents/RNAseq/read_counts_May/reps_2and3/")
#import list of library names
lib_list = list(open("/home/chuck/Documents/RNAseq/lib_names_2and3.txt", "r"))

##Initialize dataframes
initializer = pd.DataFrame(pd.read_csv("/home/chuck/Documents/RNAseq/read_counts_May/reps_2and3/LIB_187_S16_L007_May_readcounts.txt", sep="\t"))
sense_genes = pd.DataFrame(initializer["genes"])
antisense_genes = pd.DataFrame(initializer["genes"])
combined_genes = pd.DataFrame(initializer["genes"])

#Create dictionaries to rename the genes to _sense and _antisense
sense_rename = {}
antisense_rename = {}
for gene in sense_genes.iterrows():
    sense_rename[gene[1][0]] = gene[1][0] + "_sense"
    antisense_rename[gene[1][0]] = gene[1][0] + "_antisense"

#Build the sense and antisense dataframes
for lib in lib_list:
    for file in files:
        if lib.strip() in file:
            lib_name = file.split("_May")[0]
            read_counts = pd.DataFrame(pd.read_csv(counts + "/" + file, sep ="\t"))
            sense_genes[lib_name] = read_counts["sense"]
            antisense_genes[lib_name] = read_counts["antisense"]

#Set the genes as index
sense_genes = sense_genes.set_index("genes")
antisense_genes = antisense_genes.set_index("genes")

#Rename the genes based on the dictionaries from above
##Chose to set the genes as index and use df.rename becausae df.replace on a column is much slower
sense_genes = sense_genes.rename(sense_rename)
antisense_genes = antisense_genes.rename(antisense_rename)

##Create a dataframe that is a combination of alternating sense --> antisense coverage
combined = pd.DataFrame(columns = sense_genes.columns.values)
for sense, antisense in zip(sense_genes.iterrows(), antisense_genes.iterrows()):
    combined = combined.append(sense[1])
    combined = combined.append(antisense[1])

sense_genes.to_csv("/home/chuck/Documents/RNAseq/read_counts_May/sense_readcounts.txt", sep="\t", index=True)
antisense_genes.to_csv("/home/chuck/Documents/RNAseq/read_counts_May/antisense_readcounts.txt", sep="\t", index=True)
combined.to_csv("/home/chuck/Documents/RNAseq/read_counts_May/combined_readcounts.txt", sep="\t", index=True)


