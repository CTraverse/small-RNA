#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 11:02:32 2018

@author: chuck
"""

"""
This script will take each sRNA/antisense RNA and randomly shuffle its sequence to preseve the base composition
Randomly using a program that randomly shuffles the sequences but preserves dinucleotide frequency
These shuffled sequences will be used to calculate the folding energy as a baseline, random case
"""

import subprocess
from Bio import SeqIO
import os
import numpy as np

#Open up the sequences to be shuffled
sequences = list(SeqIO.parse('/home/chuck/Documents/RNAseq/master_sRNA_antisense_sequences.fasta', 'fasta'))
#Open a file that the folding energy of the shuffled sequences will be written to
RNALfold_stats = open('/home/chuck/Documents/RNAseq/RNAfolding/permutations_statistics_np.txt', 'w')

#uShuffle the sequences -- preserves dinucleotide frequencies
for sequence in sequences:
    #The program generates a new file for each shuffled sequence, so write each file to a new directory
    new_file = open('/home/chuck/Documents/RNAseq/RNAfolding/permutations_1000/' + sequence.description + "_permutations_1k.fasta", "w")
    #Build the uShuffle command that will be executed. We will generate 1000 permutations per sequences
    uShuffle_command = ['/home/chuck/Documents/RNAseq/RNAfolding/uShuffle/main.exe', '-s', str(sequence.seq), '-n', '1000', '-k', '2', '-seed', '137']
    #Run uShuffle as a subprocess
    permutations = list(subprocess.check_output(uShuffle_command).decode("utf-8").split("\n"))
    #Delete the last item in the permutations, which is just metadata we don't need
    del permutations[len(permutations) - 1]
    #Write the 1000 permutations to a new file
    for permutation in permutations:
        new_file.write(">%s_%s\n%s\n" % (sequence.description, str(permutations.index(permutation) + 1), permutation))

        
#Determine the folding energy of all of the shuffled RNA sequences. These use the files we just generated in the previous loop
for RNA_permutations in os.listdir('/home/chuck/Documents/RNAseq/RNAfolding/permutations_numpy/'):
    #Open the folding energy file that the results will be written to
    RNALfold_out = open('/home/chuck/Documents/RNAseq/RNAfolding/permutations_numpy/%s' % (RNA_permutations[:-6] + "_folding.txt"), "w")
    #Build the command for the RNA folding energy program
    RNALfold_command = ["RNALfold", '--infile=/home/chuck/Documents/RNAseq/RNAfolding/permutations_numpy/%s' % (RNA_permutations)]
    #Run RNAfold
    RNALfold = list(subprocess.check_output(RNALfold_command).decode("utf-8").split("\n"))
    
    #Loop through folding energies from each permutation
    #RNAfold outputs in a messy format, so process the format and just grab the folding energy
    #We don't care about the actual folding structure in the output
    free_energy_list = []
    for fold in RNALfold:
        if ">" in fold:
            current_permutation = fold
        if fold[0:2] == " (":
            current_energy = fold.replace("(", "").replace(" ", "").replace(")", "")
            RNALfold_out.write("%s\t%s\n" % (current_permutation, current_energy))
            free_energy_list.append(-1*float(current_energy.replace("-", "")))
            
    #Write the mean folding energies for the permutations of each original shuffled sequence
    RNALfold_stats.write("%s\t%s\t%s\t%s\t%s\n" % (RNA_permutations[:-6], np.mean(free_energy_list), np.std(free_energy_list), 
                                  np.mean(free_energy_list) + np.std(free_energy_list), 
                                  np.mean(free_energy_list) - np.std(free_energy_list)))
    
RNALfold_out.close()
