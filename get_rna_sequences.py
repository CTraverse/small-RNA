#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 15:47:07 2018

@author: chuck
"""
from Bio.Seq import Seq

"""
This script will pull all of the sequences of the antisense/sRNAs from the genome fasta file.
An extra 100 bases of sequence will be added to the region directly before the RNA starts for
downstream promoter detection.
"""

#Set up reference

reffile = "/home/chuck/Documents/RNAseq/REL606.fna"
reference = open(reffile,"r")
reference.readline()
ref = reference.readlines()
reference = "".join(ref)
reference = reference.replace("\n","") #Remove new lines and and make reference one long string
reference = " " + reference #Add throw-away character to make indexing easier

additional_region = 100

master = list(open("/home/chuck/Documents/RNAseq/master_sRNA_antisense.txt", "r"))
outfile = open("/home/chuck/Documents/RNAseq/master_sRNA_antisense_sequences_promoters.fasta", "w")

#Loop through each putative RNA in the master file and pull the sequence of where expression starts
##and ends. Also grab an extra 100 bases before the expression starts for promoter detection
for rna in master:
    rna_split = rna.strip().split("\t")
    rna_start, rna_stop = int(rna_split[2]), int(rna_split[3])
    
    #If the RNA isn't on the reverse strand, grab from the forward strand
    #If the RNA is on the reverse strand, reverse complement the forward strand into the reverse strand
    ##Then grab the sequences in the expression region and the extra 100 bases
    if 'reverse' not in rna:
        sequence = reference[int(rna_start) - additional_region:int(rna_start)]
        
    if 'reverse' in rna:
        sequence = reference[int(rna_stop):int(rna_stop) + additional_region]
        sequence = str(Seq(sequence).reverse_complement())

    header = ">%s_%s_%s_%s" % (rna_split[4], rna_split[5], rna_split[2], rna_split[3])
    outfile.write(header + "\n" + sequence + "\n")    
    
    
    
