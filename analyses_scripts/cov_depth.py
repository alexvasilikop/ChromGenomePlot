#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util import samtools_depths_util, chrom_selection
from collections import defaultdict
import matplotlib.pyplot as plt
import os
import argparse

# Plots coverage along chromosomes (all or selection). Window size is given as an argument. Input is a sorted BAM file
##############################################################################################################################
class Alignment():

	def __init__(self, filename):

		self.filename = filename
		#Dictionary of tuples: chromosome -> (position, depth)
		self.chromosomes_and_depths = defaultdict(lambda: [])
		#Dictionary of dictionaries: chromosome -> (start/end position in window -> average depth for window)
		self.chromosomes_windows_depths = defaultdict(lambda: {})
		#Dictionary with lengths of chromosomes
		self.chromosome_lengths = defaultdict(lambda: {})
		self.genome_size=0
		self.average_depth=0

	def read_depths_file(self, depths_file, selected_chrom_to_use):

		depth_total=0
		(depth_total, self.chromosomes_and_depths, self.genome_size)=samtools_depths_util.read_depths_file_util(depths_file, depth_total, selected_chrom_to_use, self.chromosomes_and_depths, self.genome_size)
		self.average_depth=depth_total/self.genome_size

	def read_alignment(self, selected_chrom_to_use):

		#If depths file has been generated in previous run read this file. Otherwise run samtools depth
		depths_file="depths.txt"

		if os.path.isfile(depths_file):
			print("Found existing depths file! SAMtools depth function will be skipped...")
			print("Using previously inferred depths.txt file...")
			self.read_depths_file(depths_file, selected_chrom_to_use)

		else:
			depths_file="depths.txt"
			samtools_depths_util.run_samtools(self.filename, depths_file)
			self.read_depths_file(depths_file, selected_chrom_to_use)

	def get_average_coverages_chromosomes(self, bin_size):

		windowstep=bin_size

		for chromosome in self.chromosomes_and_depths.keys():

			print("Working on chromosome: "+chromosome+" ...")
			start_window = 1
			depth_total_window = 0
			count_bases_window = 0
			self.chromosome_lengths[chromosome] = len(self.chromosomes_and_depths[chromosome])

			for p, d in self.chromosomes_and_depths[chromosome]:

				depth_total_window += d
				count_bases_window += 1

				if p%windowstep==0:

					average_depth = int(depth_total_window/windowstep)
					#add depth of start position
					self.chromosomes_windows_depths[chromosome][start_window] = average_depth
					#add depth of end position
					self.chromosomes_windows_depths[chromosome][p] = average_depth
					depth_total_window = 0
					start_window= p+1
					count_bases_window = 0

				elif p == self.chromosome_lengths[chromosome]:

					average_depth = int(depth_total_window/count_bases_window)
					self.chromosomes_windows_depths[chromosome][start_window] = average_depth
					self.chromosomes_windows_depths[chromosome][p] = average_depth

###############################################################################################################################
class Plot():

	def __init__(self, fig_name):
		self.fig_name=fig_name

	def plot_coverages_chromosomes(self, species, chromosomes_windows_depths, chrom_lengths, no_fill):

		plt.rcParams["font.family"]= "Arial"
		fig, axs = plt.subplots(nrows=len(chromosomes_windows_depths), ncols=1, figsize=(10,12), sharey=True, sharex=True, tight_layout=True)
		no=0

		if len(chromosomes_windows_depths.keys())>1:

			for i in chromosomes_windows_depths.keys():
				tuples = sorted(chromosomes_windows_depths[i].items()) # sorted by key, return a list of tuples
				x, y = zip(*tuples) # unpack a list of pairs into two tuples
				plt.locator_params(axis='y', nbins=4)
				axs[no].plot(x, y, color="steelblue")
				if not no_fill:
					axs[no].fill_between(x=x, y1=y, color='blue', alpha=.1)
				axs[no].locator_params(axis='y', nbins=4)
				axs[no].set_title(i, fontsize=8)
				axs[no].grid()
				no+=1

		else:
			tuples = sorted(chromosomes_windows_depths[list(chromosomes_windows_depths.keys())[0]].items()) # sorted by key
			x, y = zip(*tuples) # unpack a list of pairs into two tuples
			plt.locator_params(axis='y', nbins=4)
			axs.plot(x, y, color="steelblue")
			if not no_fill:
				axs.fill_between(x=x, y1=y, color='blue', alpha=.1)
			axs.locator_params(axis='y', nbins=4)
			axs.set_title(list(chromosomes_windows_depths.keys())[0], fontsize=8)
			axs.grid()

		interval=max(chrom_lengths.values())/8
		plt.gca().xaxis.set_major_locator(plt.MultipleLocator(interval))
		current_values = plt.gca().get_xticks().tolist()
		plt.gca().set_xticks(current_values)
		plt.gca().set_xticklabels(['{:.1f}'.format(int(x)/1000000) for x in current_values])
		plt.xlim(0,max(chrom_lengths.values())+500000)

		fig.supxlabel("Position on chromosome (Mbp)")
		fig.supylabel("Coverage depth")
		if len(chromosomes_windows_depths.keys())>1:
			fig.suptitle('Depth of coverage along chromosomes'+" - "+species, fontsize=10)
		else:
			fig.suptitle('Depth of coverage along chromosome '+list(chromosomes_windows_depths.keys())[0]+" - "+species, fontsize=10)
		fig.savefig(self.fig_name)
		plt.show()

###############################################################################################################################################################
def main():

	print("\n########################################################################################################################")
	print("Estimating coverage depth for each chromosome window and plotting coverage depth along the chromosomes of the reference genome...")
	print("##########################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot coverage depth along chromosomes from sorted BAM file. Requires SAMtools on path.")
	parser.add_argument("in_alignment", help="Genome alignment file in BAM format (sorted)")
	parser.add_argument("out_plot", help="Output plot with suffix (e.g., *.png, *.svg)")
	parser.add_argument("bin_size", help="Size of bins (bp) for coverage depth plot. Average coverage depth is plotted for each bin")
	parser.add_argument("species", help="Species name or library for the title of the plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\" )")
	parser.add_argument('--no_fill', action='store_true', help="Use this flag to indicate no color filling between the lineplot and the x axis")

	parser.usage="python3 chromgenomeplot.py cov_depth [positional_arguments]"
	args = parser.parse_args()
	selected_chrom_to_use=chrom_selection.chromosome_selection(args.sel_chrom)

	#Assign bam/sam file as an instance of the Alignment Class
	my_alignment=Alignment(filename=args.in_alignment)

	#read BAM file and calculate per base coverages -> generates depths.txt file
	print("Generating and reading depths file..\n")
	my_alignment.read_alignment(selected_chrom_to_use)

	#calculates coverage per window along the chromosomes
	print("Calculating average coverage depth per window size of "+str(args.bin_size)+" bp"+"\n")
	my_alignment.get_average_coverages_chromosomes(int(args.bin_size))

	#Print genome size (selected chtomosomes and average coverage depth)
	print("Genome size from depths file (whole-genome or selected chromosomes):\n"+str(my_alignment.genome_size)+"bp\n")
	print("Average coverage depth of selected chromosomes (whole genome or selected chromosomes):\n"+str(int(my_alignment.average_depth))+" X"+"\n")

	#plots coverages along the chromosomes
	print("Preparing the plot...\n")
	my_plot = Plot(fig_name=args.out_plot)
	my_plot.plot_coverages_chromosomes(args.species,my_alignment.chromosomes_windows_depths,my_alignment.chromosome_lengths, args.no_fill)

#############################################################################################################################################
if __name__ == '__main__':
	main()
