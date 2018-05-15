#-*- coding: utf-8 -*-
"""
Created on Mon Feb  5 14:45:37 2018

@author: chuck
"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.lines import Line2D
import os
from Bio import SeqIO

libraries = list(open("../lib_sample_names_reps_2and3.txt", "r"))

coding_genes = list(SeqIO.parse("/stor/work/Ochman/chuck/LTEE_RNAseq/Reads/REL606_coding_genes.fasta", "fasta"))

coding_dict = {}

for gene in coding_genes:
    
    if "pseudo=true" not in gene.description: #Ignore psuedogenes
        gene_split = gene.description.split("] ")
    
        if "complement" not in gene.description: #Gene is on forward strand
            forward = "Sense"
            reverse = "Antisense"
        else: #Gene is on reverse strand
            forward = "Antisense"
            reverse = "Sense"
            
        l = [gene_split.index(i) for i in gene_split if 'location=' in i][0]
        coordinates = gene_split[l].strip().strip("]").split("..")

        t = [gene_split.index(i) for i in gene_split if 'tag=' in i][0]
        locus_tag = gene_split[t].split("tag=")[1]

        ##filter non-numbers out, returns a list of lists of each individual number, so join them all together as a string, then turn it into an integer
        coords = [int("".join(list(filter(str.isdigit, coordinates[0])))), int("".join(list(filter(str.isdigit, coordinates[1]))))]
        coding_dict[locus_tag] = [coords[0], coords[1], forward, reverse]

print("Loading read information into memory...")
forward = pd.DataFrame(pd.read_csv("../intergenic_forward_1_to_4629812_norm.txt", sep="\t")) 
print("Loaded forward.")
reverse = pd.DataFrame(pd.read_csv("../intergenic_reverse_1_to_4629812_norm.txt", sep="\t"))
print("loaded reverse.")


for region in os.listdir("../RNAseq_sample_info_reps_2and3/"):
    print("Starting " + region)
    generation = region.split("_")[1].strip(".txt")
    ltee_file = list(open("../RNAseq_sample_info_reps_2and3/" + region, 'r'))
    DESeq_output = 'DESeq/combined_generation_below_vs_above_%s.txt' % (generation)
    output_dir = DESeq_output.strip(".txt") + "/"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    col = {}
    for lib in ltee_file:
        if "LIB" in lib:
                lib_split = lib.split("\t")
                if "above" in lib_split[2]:
                    col[lib_split[0]] = "gold"
                else:
                    col[lib_split[0]] = "royalblue"
        
    with open(DESeq_output, 'r') as DESeq:
        results = list(DESeq)
        del results[0]

        for result in results:
            
            if "NA" in result: # Coverage is too low at this gene to obtain a valid p-value
                pass
            
            else:
                
                result_split = result.strip().strip().split("\t")
                locus_tag = result_split[0]
                p = float(str(result_split[5]))
                
                if p <= 0.05 and "anti" in locus_tag:
                    coordinates = coding_dict[locus_tag.split("_anti")[0]]
                    intergenic_start = int(coordinates[0])
                    intergenic_stop = int(coordinates[1])
                    print("Graphing %s to %s..." % (intergenic_start, intergenic_stop))               
                    
                    width = 0.5
                    fontsize = 8
                    
                    plt.figure(1)
                    plt.subplot(211)
                    
                    region_to_plot = pd.DataFrame(forward[intergenic_start-500:intergenic_stop+500])
        
                    ax = plt.subplot(211)
                    ax.spines["top"].set_linewidth(width) 
                    ax.spines["bottom"].set_linewidth(width)   
                    ax.spines["right"].set_linewidth(width)    
                    ax.spines["left"].set_linewidth(width)  
                    ax.xaxis.set_tick_params(width=width)
                    ax.yaxis.set_tick_params(width=width)
                    
                    for lib in libraries:
                        lib = lib.strip()
                        plt.plot(region_to_plot["Locations"].astype("float"), region_to_plot[lib].astype("float"), color=col[lib], linewidth = width)
                    
                    ### SET AXIS LABELS ####
                    plt.ylabel("Standardized Expression")
                    plt.xlabel("Genomic Location (Forward Strand -- %s transcription)" % (coordinates[2]))
                    
                    ### SET VERTICAL LINES ####
                    
                    plt.axvline(x=intergenic_start, color="red", linewidth = width)
                    plt.axvline(x=intergenic_stop-1, color="red", linewidth = width)
                    
                    ### SET X AXIS TICKS, GET RID OF SCIENTIFIC NOTATION, AND SEPARATE NUMBERS WITH COMMAS ####
                    ax.xaxis.set_ticks(np.arange(intergenic_start, intergenic_stop, (intergenic_stop-intergenic_start-1)))
                    ax.get_xaxis().get_major_formatter().set_useOffset(False)
                    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
                    ax.set_xlim([intergenic_start-500, intergenic_stop+500])
                    
                    f_min, f_max = plt.gca().get_ylim()
                    
                    custom_lines = [Line2D([0], [0], color="royalblue", lw=1),
                                    Line2D([0], [0], color="gold", lw=1),
                                    Line2D([0], [0], color="red", lw=1),]
                    
                    ax.legend(custom_lines, ['Below %s generations' % (generation), 'Above %s generations' % (generation), 'Gene boundaries'], loc='upper left',prop={'size': 7})
                    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                                 ax.get_xticklabels() + ax.get_yticklabels()):
                        item.set_fontsize(fontsize)
                    
                    
                    ax = plt.subplot(212)
                    
                    region_to_plot = pd.DataFrame(reverse[intergenic_start-500:intergenic_stop+500])
        
                    ### SET LINE ####
                    ax.spines["top"].set_linewidth(width) 
                    ax.spines["bottom"].set_linewidth(width)   
                    ax.spines["right"].set_linewidth(width)    
                    ax.spines["left"].set_linewidth(width)  
                    ax.xaxis.set_tick_params(width=width)
                    ax.yaxis.set_tick_params(width=width)
                    
                    
                    for lib in libraries:
                        lib = lib.strip()
                        plt.plot(region_to_plot["Locations"].astype("float"), region_to_plot[lib].astype("float"), color=col[lib], linewidth = width)
                    
                    ### SET AXIS LABELS ####
                    plt.ylabel("Standardized Expression")
                    plt.xlabel("Genomic Location (Reverse Strand -- %s transcription)"% (coordinates[3]))
                    
                    ### SET VERTICAL LINES ####
                    plt.axvline(x=intergenic_start, color="red", linewidth = width)
                    plt.axvline(x=intergenic_stop-1, color="red", linewidth = width)
                    
                    ### SET X AXIS TICKS, GET RID OF SCIENTIFIC NOTATION, AND SEPARATE NUMBERS WITH COMMAS ####
                    ax.xaxis.set_ticks(np.arange(intergenic_start, intergenic_stop, (intergenic_stop-intergenic_start-1)))
                    ax.get_xaxis().get_major_formatter().set_useOffset(False)
                    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
                    ax.set_xlim([intergenic_start-500, intergenic_stop+500])
                    
                    r_min, r_max = plt.gca().get_ylim()
                    
                    custom_lines = [Line2D([0], [0], color="royalblue", lw=1),
                                    Line2D([0], [0], color="gold", lw=1),
                                    Line2D([0], [0], color="red", lw=1),]
                    
                    ax.legend(custom_lines, ['Below %s generations' % (generation), 'Above %s generations' % (generation), 'Gene boundaries'], loc='upper left',prop={'size': 7})
    
                    plot_min = min([f_min, r_min])
                    plot_max = max([f_max, r_max])
                    
                    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                                 ax.get_xticklabels() + ax.get_yticklabels()):
                        item.set_fontsize(fontsize)
                    
                    
                    plt.subplot(211).set_ylim([plot_min, plot_max])
                    
                    plt.subplot(212).set_ylim([plot_min, plot_max])
                    plt.tight_layout()
                    plt.savefig(output_dir + "intergenic_region_%s_to_%s_norm.pdf" % (str(intergenic_start), str(intergenic_stop)), format="pdf")
                    plt.clf()
                    plt.close()

