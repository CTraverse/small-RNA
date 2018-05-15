# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 10:42:27 2018

@author: chuck
"""

import os
import pandas as pd
from sys import argv


script, reffile, cov_files_path, out_file = argv

outfile = open(out_file, 'w')

reference = open(reffile,"r")
reference.readline()
ref = reference.readlines()
reference = "".join(ref)
reference = reference.replace("\n","")
reference = " " + reference

initialized = False



features =  pd.DataFrame(pd.read_csv("/home/chuck/Documents/RNAseq/read_counts/LIB_172_S1_L007_readcounts_2.txt",
                                                            sep="\t", names = ["features", "forward", "reverse"]))["features"]
master_forward = pd.DataFrame([[features[0]+"_forward"], [features[0]+"_reverse"]], columns = ['features'])
master_reverse = pd.DataFrame([[features[0]+"_forward"], [features[0]+"_reverse"]], columns = ['features'])

del features[0]

for row in features:

    add_next = pd.DataFrame([[row+"_forward"], [row+"_reverse"]], columns = ['features'])

    master_read_counts = master_read_counts.append(add_next, ignore_index=True)


for library in os.listdir(cov_files_path):
    if initialized == False:
        current_coverage_file =  pd.DataFrame(pd.read_csv(cov_files_path + "/" + library,
                                                            sep="\t", names = ["position", "base", "foward", "reverse"]))
        
        master_forward = current_coverage_file.drop('reverse').rename(columns={'forward': library})
        
        master_reverse = current_coverage_file.drop('forward').rename(columns={'reverse': library})
        initialized = True
    else:
        master_forward
    
    
    current_library = pd.DataFrame(pd.read_csv(os.path.join("/home/chuck/Documents/RNAseq/read_counts", library), sep="\t", names = ["features", "forward", "reverse"]))
    current_library.sort_index()
    F, R = current_library.loc[0]["forward"], current_library.loc[0]["reverse"]
    
    add_next = pd.DataFrame([[F], [R]], columns = [library.split("_readcounts")[0]])

    current_library = current_library.iloc[1:]

    for index, row in current_library.iterrows():
        F, R = row["forward"], row["reverse"]

        current_feature = pd.DataFrame([[F], [R]], columns = [library.split("_readcounts")[0]])
        add_next = add_next.append(current_feature, ignore_index=True)
    new_master = [master_read_counts, add_next]
    master_read_counts[library.split("_readcounts")[0]] = add_next
    print "Done with " + library.split("_readcounts")[0]

master_read_counts.to_csv("read_counts_for_DEseq.txt", sep="\t")












