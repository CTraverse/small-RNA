#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 15:47:07 2018

@author: chuck
"""
from Bio.Seq import Seq
#from Bio.SeqUtils import GC
##This script will calculate all of the sequences of the antisense/sRNAs

reffile = "/home/chuck/Documents/RNAseq/REL606.fna"
#Set up reference
reference = open(reffile,"r")
reference.readline()
ref = reference.readlines()
reference = "".join(ref)
reference = reference.replace("\n","") #Remove new lines and and make reference one long string
reference = " " + reference #Add throw-away character to make indexing easier

additional_region = 100

master = list(open("/home/chuck/Documents/RNAseq/master_sRNA_antisense.txt", "r"))
outfile = open("/home/chuck/Documents/RNAseq/master_sRNA_antisense_sequences_promoters.fasta", "w")

for rna in master:
    rna_split = rna.strip().split("\t")
    rna_start, rna_stop = int(rna_split[2]), int(rna_split[3])
    if 'reverse' not in rna:
        
        sequence = reference[int(rna_start) - additional_region:int(rna_start)]
    if 'reverse' in rna:
        sequence = reference[int(rna_stop):int(rna_stop) + additional_region]

        sequence = str(Seq(sequence).reverse_complement())

    header = ">%s_%s_%s_%s" % (rna_split[4], rna_split[5], rna_split[2], rna_split[3])
    outfile.write(header + "\n" + sequence + "\n")
#    print(header, len(sequence), GC(Seq(sequence)))
    
    
    
    