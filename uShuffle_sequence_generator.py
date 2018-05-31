#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 11:02:32 2018

@author: chuck
"""

import subprocess
from Bio import SeqIO
import os
import numpy as np

sequences = list(SeqIO.parse('/home/chuck/Documents/RNAseq/master_sRNA_antisense_sequences.fasta', 'fasta'))

RNALfold_stats = open('/home/chuck/Documents/RNAseq/RNAfolding/permutations_statistics_np.txt', 'w')

preserve_dinuc = 0

### This script will take each sRNA/antisense RNA and randomly shuffle its sequence to preseve the base composition
#### Two options: randomly shuffle, or use a more sophisticated program that randomly shuffles but preserves dinucleotide frequency

if preserve_dinuc == 0: ###numpy shuffler 

    for sequence in sequences:
        seq = list(str(sequence.seq))
        new_file = open('/home/chuck/Documents/RNAseq/RNAfolding/permutations_numpy/' + sequence.description + "_permutations_np.fasta", "w")
    
        for i in range(1000):
            new_seq = "".join(np.random.RandomState(seed=i).permutation(seq))
            new_file.write(">%s_%s\n%s\n" % (sequence.description, i + 1, new_seq))
            
            
else: #### uShuffle method -- preserves dinucleotide frequencies
        
    for sequence in sequences:
        new_file = open('/home/chuck/Documents/RNAseq/RNAfolding/permutations_1000/' + sequence.description + "_permutations_1k.fasta", "w")
        
        uShuffle_command = ['/home/chuck/Documents/RNAseq/RNAfolding/uShuffle/main.exe', '-s', str(sequence.seq), '-n', '1000', '-k', '2', '-seed', '137']
    
        permutations = list(subprocess.check_output(uShuffle_command).decode("utf-8").split("\n"))
        del permutations[len(permutations) - 1]
        for permutation in permutations:
            new_file.write(">%s_%s\n%s\n" % (sequence.description, str(permutations.index(permutation) + 1), permutation))

        
### Determine the folding energy of all of the shuffled RNA sequences
for RNA_permutations in os.listdir('/home/chuck/Documents/RNAseq/RNAfolding/permutations_numpy/'):
    RNALfold_out = open('/home/chuck/Documents/RNAseq/RNAfolding/permutations_numpy/%s' % (RNA_permutations[:-6] + "_folding.txt"), "w")
    RNALfold_command = ["RNALfold", '--infile=/home/chuck/Documents/RNAseq/RNAfolding/permutations_numpy/%s' % (RNA_permutations)]
    RNALfold = list(subprocess.check_output(RNALfold_command).decode("utf-8").split("\n"))
    
    free_energy_list = []
    for fold in RNALfold:
        if ">" in fold:
            current_permutation = fold
        if fold[0:2] == " (":
            current_energy = fold.replace("(", "").replace(" ", "").replace(")", "")
            RNALfold_out.write("%s\t%s\n" % (current_permutation, current_energy))
            free_energy_list.append(-1*float(current_energy.replace("-", "")))
    RNALfold_stats.write("%s\t%s\t%s\t%s\t%s\n" % (RNA_permutations[:-6], np.mean(free_energy_list), np.std(free_energy_list), 
                                  np.mean(free_energy_list) + np.std(free_energy_list), 
                                  np.mean(free_energy_list) - np.std(free_energy_list)))
    
#    print("%s\t%s\t%s\t%s\t%s" % (RNA_permutations[:-6], np.mean(free_energy_list), np.std(free_energy_list), 
#                                    np.mean(free_energy_list) + np.std(free_energy_list), 
#                                    np.mean(free_energy_list) - np.std(free_energy_list)))

RNALfold_out.close()