# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 16:58:56 2017

@author: chuck
"""
import gzip
import numpy as np
from sys import argv
import pandas as pd
import time
from Bio import SeqIO

script, samfile_path, phred_cutoff, reffile, genes_file, out_file = argv

"""
This script builds a custom RNAseq read merger and then then determines the overall 
per-base coverage and counts the total number of reads per gene across the genome.
"""



###########################
#########Functions#########
###########################

##Find the highest quality read for each overlapping base in mate pairs
def highest_quality(mate_1_overlap_seq, mate_2_overlap_seq, mate_1_overlap_quality, mate_2_overlap_quality):
        mate_overlap_consensus = ""
        score_overlap_consensus = ""
        #Read in base calls and their quality scores
        for base_1, base_2, score_1, score_2 in zip(mate_1_overlap_seq, mate_2_overlap_seq, mate_1_overlap_quality, mate_2_overlap_quality):
            score_1_prob, score_2_prob = (1 - 10**(float(ord(score_1)-33)/-10)), (1 - 10**(float(ord(score_2)-33)/-10)) #Transform quality scores into probability scores
            
            #Combine the two reads where they overlap
            if score_1_prob >= score_2_prob:
                mate_overlap_consensus = mate_overlap_consensus + base_1
                score_overlap_consensus = score_overlap_consensus + score_1
    
            else:
                mate_overlap_consensus = mate_overlap_consensus + base_2
                score_overlap_consensus = score_overlap_consensus + score_2
                
        return mate_overlap_consensus, score_overlap_consensus


# This function is a custom read merger. The comments with dashes below are the types of reads the code is meant to resolve
def consensus(mate_1_position, mate_2_position, mate_1_sequence, mate_2_sequence, mate_1_quality, mate_2_quality):
    mate_offset = mate_2_position - mate_1_position
    mate_1_len, mate_2_len = len(mate_1_sequence), len(mate_2_sequence)
    
    if mate_offset == 0: #Reads start at same position
        
        if mate_1_len == mate_2_len: #Perfect overlap
            mate_consensus, quality_consensus = highest_quality(mate_1_sequence, mate_2_sequence, mate_1_quality, mate_2_quality)
            ########--------------------########
            ########--------------------########
                        
        elif mate_1_len > mate_2_len: #Mate 1 is longer than mate 2
            mate_consensus, quality_consensus = highest_quality(mate_1_sequence[0:mate_2_len], mate_2_sequence, mate_1_quality[0:mate_2_len], mate_2_quality)
            mate_consensus = mate_consensus + mate_1_sequence[mate_2_len:mate_1_len]
            quality_consensus = quality_consensus + mate_1_quality[mate_2_len:mate_1_len]
            ########-------------------------###
            ########--------------------########
                        
        else: # mate_2 is longer than mate 1
            mate_consensus, quality_consensus = highest_quality(mate_1_sequence, mate_2_sequence[0:mate_1_len], mate_1_quality, mate_2_quality[0:mate_1_len])
            mate_consensus = mate_consensus + mate_2_sequence[mate_1_len:mate_2_len]
            quality_consensus = quality_consensus + mate_2_quality[mate_1_len:mate_2_len]
            ########--------------------########
            ########-------------------------###
                   
            
    elif mate_offset == mate_1_len: #Mate 1 and mate 2 are exactly next to one another
        mate_consensus = mate_1_sequence + mate_2_sequence
        quality_consensus = mate_1_quality + mate_2_quality
            ########--------------------############################
            ############################--------------------########
        
        
    elif mate_offset > 0: #Mate 1 mapped before mate 2
        
        if mate_1_len == (mate_offset + mate_2_len): #They end at the same position
            mate_consensus, quality_consensus = highest_quality(mate_1_sequence[mate_offset:mate_1_len], mate_2_sequence, mate_1_quality[mate_offset:mate_1_len], mate_2_quality)
            mate_consensus = mate_1_sequence[0:mate_offset] + mate_consensus
            quality_consensus = mate_1_quality[0:mate_offset] + quality_consensus
            ####------------------------########
            ########--------------------########
            
        elif mate_1_len > (mate_offset + mate_2_len):# Mate 2 starts and ends within mate 1; just use mate 1
            mate_consensus, quality_consensus = highest_quality(mate_1_sequence[mate_offset:(mate_offset+mate_2_len)], mate_2_sequence, mate_1_quality[mate_offset:(mate_offset+mate_2_len)], mate_2_quality)
            mate_consensus = mate_1_sequence[0:mate_offset] + mate_consensus + mate_1_sequence[(mate_offset+mate_2_len):mate_1_len]
            quality_consensus = mate_1_quality[0:mate_offset] + quality_consensus + mate_1_quality[(mate_offset+mate_2_len):mate_1_len]
            ####------------------------########
            ########----------------############      
                    
        else: #They start and end in different positions
            mate_consensus, quality_consensus = highest_quality(mate_1_sequence[mate_offset:mate_1_len], mate_2_sequence[0:(mate_1_len - mate_offset)], mate_1_quality[mate_offset:mate_1_len], mate_2_quality[0:(mate_1_len - mate_offset)])
            mate_consensus = mate_1_sequence[0:mate_offset] + mate_consensus + mate_2_sequence[(mate_1_len - mate_offset):mate_2_len]
            quality_consensus = mate_1_quality[0:mate_offset] + quality_consensus + mate_2_quality[(mate_1_len - mate_offset):mate_2_len]
            ####------------------------########
            ########------------------------####     

        
    elif mate_offset < 0: #Mate 2 mapped before mate 1
        mate_offset = abs(mate_offset)
        if mate_2_len == (mate_offset + mate_1_len): #They end at the same position
            mate_consensus, quality_consensus = highest_quality(mate_1_sequence, mate_2_sequence[mate_offset:mate_2_len], mate_1_quality, mate_2_quality[mate_offset:mate_2_len])
            mate_consensus = mate_2_sequence[0:mate_offset] + mate_consensus
            quality_consensus = mate_2_quality[0:mate_offset] + quality_consensus
            ########--------------------########
            ####------------------------########  
            
        elif mate_2_len > (mate_offset + mate_1_len):# Mate 1 starts and ends within mate 2; just use mate 2
            mate_consensus, quality_consensus = highest_quality(mate_1_sequence, mate_2_sequence[mate_offset:(mate_offset+mate_1_len)], mate_1_quality, mate_2_quality[mate_offset:(mate_offset+mate_1_len)])
            mate_consensus = mate_2_sequence[0:mate_offset] + mate_consensus + mate_2_sequence[(mate_offset+mate_1_len):mate_2_len]
            quality_consensus = mate_2_quality[0:mate_offset] + quality_consensus + mate_2_quality[(mate_offset+mate_1_len):mate_2_len]
            ########----------------############                         
            ####------------------------########

        else:
            mate_consensus, quality_consensus = highest_quality(mate_1_sequence[0:(mate_2_len - mate_offset)], mate_2_sequence[mate_offset:mate_2_len], mate_1_quality[0:(mate_2_len - mate_offset)], mate_2_quality[mate_offset:mate_2_len])
            mate_consensus = mate_2_sequence[0:mate_offset] + mate_consensus + mate_1_sequence[(mate_2_len - mate_offset):mate_1_len]
            quality_consensus = mate_2_quality[0:mate_offset] + quality_consensus + mate_1_quality[(mate_2_len - mate_offset):mate_1_len]
            ########------------------------####     
            ####------------------------########


    return mate_consensus, quality_consensus


#Set up reference genome
reference = open(reffile,"r")
reference.readline()
ref = reference.readlines()
reference = "".join(ref)
reference = reference.replace("\n","") #Remove new lines and and make reference one long string
reference = " " + reference #Add throw-away character to make indexing easier

coding_sequences = list(SeqIO.parse(genes_file, 'fasta'))
read_counts = np.zeros((len(coding_sequences), 2), dtype=int)
genes_list = []
genes_info = []
#Create a list of lists for each gene. The list for each gene contains the locus tag, strand information, and genomic coordinates
for gene in coding_sequences:
    
    if "pseudo=true" not in gene.description: #Ignore psuedogenes
        gene_split = gene.description.split("] ")
    
        if "complement" not in gene.description: #Gene is on forward strand
            strand = 0
        else: #Gene is on reverse strand
            strand = 1
    
        l = [gene_split.index(i) for i in gene_split if 'location=' in i][0]
        coordinates = gene_split[l].strip().strip("]").split("..")

        t = [gene_split.index(i) for i in gene_split if 'tag=' in i][0]
        locus_tag = gene_split[t].split("tag=")[1]

        ##filter non-numbers out, returns a list of lists of each individual number, so join them all together as a string, then turn it into an integer
        coords = [int("".join(list(filter(str.isdigit, coordinates[0])))), int("".join(list(filter(str.isdigit, coordinates[1]))))]
        genes_list.append(locus_tag)
        genes_info.append([strand, coords[0], coords[1]])
        
#Sam file flags
forward_flags = [97, 99, 145, 147]
reverse_flags = [81, 83, 161, 163]

mappings = np.zeros((len(reference), 2), dtype=int)

#Open up the alignment file (samfile) and iterate through each alignment
##Figure out if it's a good alignment, if so move forward,
##Then finds how much overlap the two reads have,
##
with gzip.open(samfile_path, 'rb') as samfile:
    current_pair = ""
    start_time = time.time()

    for line in samfile:
        line = line.decode('utf8')
        linesplit = line.split("\t")
        pair_name = linesplit[0]

        if pair_name != current_pair: #Finding paired reads
            mate_1 = linesplit
            mate_1_unsplit = line
            current_pair = pair_name
            
        elif pair_name == current_pair: #Finding paired reads
            mate_2 = linesplit
            mate_2_unsplit = line
            CIGAR_check = mate_1[5] + mate_2[5]
            
            if "I" not in CIGAR_check and "D" not in CIGAR_check and "S" not in CIGAR_check: #Only use good read alignments        
                flags = [int(mate_1[1]), int(mate_2[1])]
                mate_1_pos, mate_2_pos = int(mate_1[3]), int(mate_2[3])
                mate_1_seq, mate_2_seq = mate_1[9], mate_2[9]
                mate_1_Q, mate_2_Q = mate_1[10], mate_2[10]
                
                if mate_1_pos <= mate_2_pos:
                    first_pos, second_pos = mate_1_pos, mate_2_pos
                    first_len= len(mate_1_seq)
                    
                else:
                    first_pos, second_pos = mate_2_pos, mate_1_pos
                    first_len = len(mate_2_seq)
                    
                if (second_pos - first_pos) <= first_len: #Only accept reads that have some degree of overlap or that mapped directly next to one another
                    read_start = min([mate_1_pos, mate_2_pos])
                    
                    if flags[0] in forward_flags and flags[1] in forward_flags or flags[0] in reverse_flags and flags[1] in reverse_flags: #Perfect mapping to forward strand
                       
                        if flags[0] in forward_flags: #Which strand is the read mapping to?
                            read_strand = 0
                        else:
                            read_strand = 1
                        
                        #Get consensus of reads
                        mate_consensus, quality_consensus = consensus(mate_1_pos, mate_2_pos, mate_1_seq, mate_2_seq, mate_1_Q, mate_2_Q)
                        read_end = read_start + len(mate_consensus)
                        i = 0
                        below_threshold = 0
                        above_threshold = 0
                        
                        #Determine the overall quality of the overlapped reads based on quality scores
                        while i < len(mate_consensus):
                            phred = -1*(10*(float(ord(quality_consensus[i])-33)/-10)) #reduced this equation to avoid unnecessary calculations
                            if float(phred_cutoff) <= phred:
                                above_threshold +=1
                            else:
                                below_threshold += 1
                            i += 1
                        
                        #Only keep the reads that are at least 75% high quality
                        proportion_above_threshold = above_threshold / (above_threshold + below_threshold)
                        
                        if proportion_above_threshold >= 0.75:
                            mappings[read_start:read_end,  read_strand] += 1
                            
                            #If the read passed the quality threshold, count its overall coverage in its location in the genome
                            gene_index = 0 #Index counter
                            for gene in genes_info:
                                gene_strand, start, stop = gene
                                
                                if read_start >= start and read_start <= stop or read_end >= start and read_end <= stop: #Check if at least part of the read maps to the coding region
                                        
                                    if gene_strand == 0 and  read_strand == 0 or gene_strand == 1 and  read_strand == 1: #Sense transcription
                                        read_counts[gene_index, 0] += 1
                                    elif gene_strand == 0 and  read_strand == 1 or gene_strand == 1 and  read_strand == 0: #Antisense transcription
                                        read_counts[gene_index, 1] += 1
                                    break
                                
                                else:
                                    gene_index += 1
                                        
                    else: # The mapping is odd, so we won't use this pair of reads
                        pass
                    
#Turn the coverage and counts data into a pandas dataframe                 
mappings = pd.DataFrame(mappings, columns=["forward", "reverse"])
read_counts = pd.DataFrame(read_counts, columns=["sense", "antisense"])

#Build up the whole read counts and then write it to a file
master_read_counts = pd.DataFrame({"genes": genes_list})
master_read_counts["sense"] = read_counts["sense"]
master_read_counts["antisense"] = read_counts["antisense"]
master_read_counts.to_csv(out_file + "_readcounts.txt", sep="\t", index=False)

#Build up the whole coverage dataframe and then write it to a file
master_dataframe = pd.DataFrame(np.array(range(0,len(reference))), columns=["position"])
master_dataframe["base"] = pd.Series(list(reference)).astype("str")
master_dataframe["forward"] = mappings["forward"]
master_dataframe["reverse"] = mappings["reverse"]
master_dataframe = master_dataframe.iloc[1:] #Drop the throw-away row
master_dataframe.to_csv(out_file + "_coverage.txt", sep="\t", index=False)

