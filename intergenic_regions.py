# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 15:26:45 2017

@author: chuck
"""
from Bio import SeqIO
genome = list(SeqIO.parse('/home/chuck/Documents/RNAseq/REL606_2.gff', 'fasta'))

intergenic_write = open('/home/chuck/Documents/RNAseq/REL606_genic_and_intergenic_3.txt', 'w')

previous_end = genome[0].description.split('location=')[1].strip(']\n').split('..')[1]

intergenic_list = []
intergenic_dict = {}
del genome[0]
i = 0
intergenic_count = 0
gene_count = 0
for gene in genome:
    if 'complement' not in gene.description:
        locations = gene.description.split('location=')[1].strip(']\n')
        
    else:
        locations = gene.description.split('complement(')[1].strip(')]\n')

    gene_count += 1
    current_start = locations.split('..')[0].strip('<').strip('>')
    dist = int(current_start)-int(previous_end)    
    intergenic_list.append(int(current_start)-int(previous_end))
    if dist not in intergenic_dict:
        intergenic_dict[dist] = 1
    else:
        intergenic_dict[dist] += 1
        
    if int(current_start) - int(previous_end) > 100:
        intergenic_count += 1
        i += 1
        intergenic_write.write("%s\t%s\t%s\n" % ("intergenic_" + str(intergenic_count), str(int(previous_end)+25), str(int(current_start)-25)))
    intergenic_write.write("%s\t%s\t%s\n" % ("gene_" + str(gene_count), current_start, locations.split('..')[1].strip('<').strip('>')))

    previous_end = locations.split('..')[1].strip('<').strip('>')

new_list = sorted(list(set(intergenic_list)))

for item in new_list:
    print item, intergenic_dict[item]
print i

intergenic_write.close()