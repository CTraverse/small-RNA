# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 11:11:54 2017

@author: chuck
"""

import gzip, re
from sys import argv
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


script, in_sam, out_fastq = argv

#The RNAseq dataset used in this study was contaminated with fish DNA during the sequencing run
##This script removes all traces of the fish DNA and only keeps reads that properly aligned to the E. coli genome
##It also keeps portions of reads that aligned to the genome -- There were many fish-E. coli chimeric reads

new_fastq_1 = open(out_fastq + "clipped_1.fastq", "w")
new_fastq_2 = open(out_fastq + "clipped_2.fastq", "w")

#Below are sam file flags that signify good alignments
forward_flags = [97, 99, 145, 147]
reverse_flags = [81, 83, 161, 163]

good_flags = forward_flags + reverse_flags

with gzip.open(in_sam, 'rb') as samfile:
    current_pair = ""
    
    for line in samfile:
        linesplit = line.split("\t")
        check_flag = int(linesplit[1])
        
        if check_flag in good_flags: #Make sure the reads aligned properly
            pair_name = linesplit[0]
            
            if pair_name != current_pair: #The current read is not paired with the previous read
                mate_1 = linesplit
                current_pair = pair_name
                
            elif pair_name == current_pair:#The current read is paired with the previous read
                mate_2 = linesplit
                CIGAR_check = mate_1[5] + mate_2[5]
                current_pair = ""
                
                if "I" in CIGAR_check or "D" in CIGAR_check: #Get rid of improperly aligning reads
                    pass
                else:

                    old_1_seq, old_2_seq = mate_1[9], mate_2[9]
                    old_1_qual, old_2_qual = mate_1[10], mate_2[10]
            
                    CIGAR_1_split = re.split('(\d+)', mate_1[5])
                    CIGAR_2_split = re.split('(\d+)', mate_2[5])
                    del CIGAR_1_split[0]
                    del CIGAR_2_split[0]
                    
                    new_1_seq = old_1_seq
                    new_1_qual = old_1_qual
                    new_2_seq = old_2_seq
                    new_2_qual = old_2_qual 
                    
                    ###Clip off the parts of the reads that didn't align to the genome###
                    if CIGAR_1_split[1] == "S": 
                        bad_sequence = int(CIGAR_1_split[0]) #Sequence that didn't align
                        new_1_seq = new_1_seq[bad_sequence:len(new_1_seq)] #Only keep the sequence that aligned
                        new_1_qual = new_1_qual[bad_sequence:len(new_1_qual)] #Only keep the Qscores that aligned
        
                    if CIGAR_1_split[len(CIGAR_1_split) - 1] == "S":
                        bad_sequence = int(CIGAR_1_split[len(CIGAR_1_split)-2]) #Sequence that didn't align
                        new_1_seq = new_1_seq[0:len(new_1_seq)-bad_sequence] #Only keep the sequence that aligned
                        new_1_qual = new_1_qual[0:len(new_1_qual)-bad_sequence] #Only keep the Qscores that aligned
                        
                        
                    if CIGAR_2_split[1] == "S":
                        bad_sequence = int(CIGAR_2_split[0])#Sequence that didn't align
                        new_2_seq = new_2_seq[bad_sequence:len(new_2_seq)]#Only keep the sequence that aligned
                        new_2_qual = new_2_qual[bad_sequence:len(new_2_qual)] #Only keep the Qscores that aligned
        
                    if CIGAR_2_split[len(CIGAR_2_split) - 1] == "S":
                        bad_sequence = int(CIGAR_2_split[len(CIGAR_2_split)-2]) #Sequence that didn't align
                        new_2_seq = new_2_seq[0:len(new_2_seq)-bad_sequence] #Only keep the sequence that aligned
                        new_2_qual = new_2_qual[0:len(new_2_qual)-bad_sequence] #Only keep the Qscores that aligned

                    
                    flag_1, flag_2 = int(mate_1[1]), int(mate_2[1])
                    
                    #Write the good reads to new files. The code below converts the sam file into a fastq file
                    if flag_1 in forward_flags and flag_2 in forward_flags: #Write a new fastq file. Make sure that the forward/reverse reads are in the correct orientation
                        new_2_seq = Seq(new_2_seq, IUPAC.unambiguous_dna).reverse_complement()
                        new_fastq_1.write("@%s\n%s\n+\n%s\n" % (linesplit[0], new_1_seq, new_1_qual))
                        new_fastq_2.write("@%s\n%s\n+\n%s\n" % (linesplit[0], new_2_seq, new_2_qual)) 
                        
                    elif flag_1 in reverse_flags and flag_2 in reverse_flags:
                        new_1_seq = Seq(new_1_seq, IUPAC.unambiguous_dna).reverse_complement()
                        new_fastq_1.write("@%s\n%s\n+\n%s\n" % (linesplit[0], new_1_seq, new_1_qual))
                        new_fastq_2.write("@%s\n%s\n+\n%s\n" % (linesplit[0], new_2_seq, new_2_qual))  
                  
                    else:
                        pass
new_fastq_1.close()
new_fastq_2.close()





















