#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util import samtools_depths_util, chrom_selection
from cov_depth import Alignment
from feature import Assembly_FASTA
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import argparse

# Plots coverage along chromosomes (all or selection). Window size is given as an argument. Input is a provided sorted BAM file
##############################################################################################################################
class Alignment_new(Alignment):

	def __init__(self, filename):
		Alignment.__init__(self, filename)
		self.names= []
		self.average_depths = []
		self.lengths = []
		self.GC_contents = []

	def get_average_coverages_chromosomes(self, scaffolds_seqs):

		for chromosome in self.chromosomes_and_depths.keys():
			print("Working on chromosome: "+chromosome+" ...")
			depth_total = 0
			average_depth = 0
			self.chromosome_lengths[chromosome] = len(self.chromosomes_and_depths[chromosome])

			A_content=scaffolds_seqs[chromosome].upper().count("A")
			T_content=scaffolds_seqs[chromosome].upper().count("T")
			G_content=scaffolds_seqs[chromosome].upper().count("G")
			C_content=scaffolds_seqs[chromosome].upper().count("C")

			for p, d in self.chromosomes_and_depths[chromosome]:
				depth_total += d

				if p == self.chromosome_lengths[chromosome]:
					average_depth = int(depth_total/len(self.chromosomes_and_depths[chromosome]))
					self.average_depths.append(average_depth)
					self.GC_contents.append(((G_content+C_content)/(G_content+C_content+A_content+T_content))*100)
					print(f"GC content: {((G_content+C_content)/(G_content+C_content+A_content+T_content))*100:.2f} %")
					self.names.append(chromosome)
					length='{:.0f}'.format(len(str(scaffolds_seqs[chromosome]).upper().replace("N", "")))
					self.lengths.append(length)

	def prepare_data_structure_for_plot(self):
		gc_dataframe = pd.DataFrame({'Name':                    self.names, 
			                      'Average coverage depth':     self.average_depths, 
			                      'Length (Mb)':                self.lengths, 
			                      'GC content (%)':             self.GC_contents})

		#interpret length as float
		gc_dataframe['Length (Mb)'] = gc_dataframe['Length (Mb)'].astype(float)
		return gc_dataframe

###############################################################################################################################
class Plot():

	def __init__(self, fig_name):
		self.fig_name=fig_name

	def plot_coverages_vs_GC(self, df):

		plt.rcParams["font.family"]= "Arial"

		if df["Name"].count() <= 20:
			fig=plt.figure(figsize=(12,12), tight_layout=True)
			sns.set_style("darkgrid")
			sns.scatterplot(data=df, x="GC content (%)", y="Average coverage depth", size="Length (Mb)", hue="Name", alpha=.6)
			lgd = plt.legend(loc="upper right",bbox_to_anchor=(1.15,1.00), title="Length (Mb)")
			fig.savefig(self.fig_name, bbox_extra_artists=(lgd))
			plt.show()
		else:
			fig=plt.figure(figsize=(12,12), tight_layout=True)
			sns.set_style("darkgrid")
			sns.scatterplot(data=df, x="GC content (%)", y="Average coverage depth", size="Length (Mb)", alpha=.6)
			lgd = plt.legend(loc="upper right", bbox_to_anchor=(1.15,1.00), title="Length (Mb)")
			fig.savefig(self.fig_name, bbox_extra_artists=(lgd))
			plt.show()

###############################################################################################################################################################
def main():

	print("\n########################################################################################################################")
	print("Plotting average coverage depth vs GC content of chromosomes of the reference genome...")
	print("##########################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot coverage depth vs GC content and size of chromosomes in the genome assembly")
	parser.add_argument("in_alignment", help="Genome reads' alignment file in BAM format (sorted)")
	parser.add_argument("in_assembly", help="Genome assembly file used for the mapping")
	parser.add_argument("out_plot", help="Output plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\" )")

	parser.usage="python3 chromgenomeplot.py cov_vs_GC [positional_arguments]"
	args = parser.parse_args()
	selected_chrom_to_use=chrom_selection.chromosome_selection(args.sel_chrom)

	#Assign bam/sam file as an instance of the Alignment Class
	my_alignment=Alignment_new(filename=args.in_alignment)

	#read BAM file and calculate per base coverages -> generates depths.txt file
	print("Generating and reading Coverage depths file..\n")
	my_alignment.read_alignment(selected_chrom_to_use)
	print("Finished processing depths file..\n")

	#Parsing assembly fasta
	print("Processing Assembly fasta file...")
	my_assembly = Assembly_FASTA(filename=args.in_assembly)
	print("Reading Assembly fasta file...\n")
	my_assembly.read_fasta(selected_chrom_to_use)
	print("Finished reading Assembly fasta file...\n"+200*"-"+"\n")
	#calculates coverage per window along the chromosomes
	print("Calculating average coverage, GC content and size of chromosomes...\n")
	my_alignment.get_average_coverages_chromosomes(my_assembly.scaffolds_seqs)

	#Print genome size (selected chtomosomes and average coverage depth)
	print("Genome size from depths file (whole genome or selected chromosomes):\n"+str(my_alignment.genome_size)+" bp\n")
	print("Average coverage depth (whole genome):\n"+str(int(my_alignment.average_depth))+" X"+"\n")

	#plots coverages along the chromosomes
	print("Preparing the plot...\n")
	my_plot = Plot(fig_name=args.out_plot)
	my_plot.plot_coverages_vs_GC(df=my_alignment.prepare_data_structure_for_plot())

#############################################################################################################################################
if __name__ == '__main__':
	main()
