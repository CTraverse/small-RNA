# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 15:26:45 2017

@author: chuck
"""


"""
This script will loop through all genes in the genome and find the regions between the genes, termed 'intergenic regions'.
It builds in an extra 100 bases on either end of the intergenic region for analysis of the surrounding sequence
"""

from Bio import SeqIO

#import the genome file with gene coordinates
genome = list(SeqIO.parse('/home/chuck/Documents/RNAseq/REL606_2.gff', 'fasta'))

#open the file that the intergenic regions will be written to
intergenic_write = open('/home/chuck/Documents/RNAseq/REL606_genic_and_intergenic_3.txt', 'w')

#Initialize previous end because there will not have been a previous gene at the start of the loop
previous_end = genome[0].description.split('location=')[1].strip(']\n').split('..')[1]

#Delete the first gene in the genome since previous end is already initialized
##Start the actual loop on the second gene
del genome[0]
intergenic_count = 0
gene_count = 0

#Loop through each gene in the genome and find the regions between each consecutive gene
for gene in genome:
    if 'complement' not in gene.description: #Determine if the gene is on the forward or reverse strand
        locations = gene.description.split('location=')[1].strip(']\n')
        
    else:
        locations = gene.description.split('complement(')[1].strip(')]\n')
    
    #Counter for which gene we're on
    gene_count += 1
    #What is the coordinate designating the start of the current gene in the loop?
    current_start = locations.split('..')[0].strip('<').strip('>')

    #Only take the intergenic regions that are greater than 100 bases in length
    #Then write the locations if each gene and each intergenic region in order
    if int(current_start) - int(previous_end) > 100:
        intergenic_count += 1
        intergenic_write.write("%s\t%s\t%s\n" % ("intergenic_" + str(intergenic_count), str(int(previous_end)+25), str(int(current_start)-25)))
        
    intergenic_write.write("%s\t%s\t%s\n" % ("gene_" + str(gene_count), current_start, locations.split('..')[1].strip('<').strip('>')))
    
    #update the previous end to be used in the next iteration
    previous_end = locations.split('..')[1].strip('<').strip('>')


intergenic_write.close()
