#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util import chrom_selection, feature_utils
from collections import defaultdict
import matplotlib.pyplot as plt
import argparse
import re

###########################################################################################################################################################################################################################
# Plots two features along chromosomes (e.g. gene content vs repetitive content) using bed files of annotations (each bed annotation must contain sorted annotations by genome position annotations)
##############################################################################################################################################################################################################

class BEDfile():

	def __init__(self, filename):

		self.filename = filename
		self.chrom_features = defaultdict(lambda: [])

	def read_bed(self, selected_chrom_to_use):

		self.chrom_features=feature_utils.bed_read_util(self.filename, self.chrom_features, selected_chrom_to_use)

########################################################################################################################################
class Assembly_FASTA():

	def __init__(self, filename):

		self.filename = filename
		self.scaffolds_seqs = defaultdict(lambda: "")
		self.chromosomes_windows = defaultdict(lambda: {})

	def read_fasta(self, selected_chrom_to_use):

		self.scaffolds_seqs=feature_utils.fasta_read_util(self.filename, self.scaffolds_seqs, selected_chrom_to_use)

	def calculate_average_bed_values_chromosomes(self, bin_size, feature_1):

		'''Provide as input window size, feature_1 and feature_2 position dictionaries {chromX: [p1, p2, ..], chromY: [p1, p2, ..]}'''
		for chromosome in self.scaffolds_seqs.keys():

			print("Working on chromosome: "+chromosome+" ...")
			#Coverage-specific values
			start_window = 1
			count_bases_window = 0
			length_chromosome = len(self.scaffolds_seqs[chromosome])
			position=0
			#account for non-overlapping annotations (return unique list of covered positions)
			unique_feature_1=list(set(feature_1[chromosome]))

			for nuc in self.scaffolds_seqs[chromosome]:
				count_bases_window += 1

				if (position+1)%bin_size==0:

					#Infer percent coverage for window
					feature_1_window  = feature_utils.total_window(start_window, position, chromosome, unique_feature_1)	
					feature_1_coverage = float(feature_1_window/bin_size)*100

					#Add depth, snp percent of start position
					self.chromosomes_windows[chromosome][start_window] = (feature_1_coverage)
					#Add depth , snp percent  of end position
					self.chromosomes_windows[chromosome][position] = (feature_1_coverage)
					start_window= position+1
					count_bases_window = 0

				#if reaching the end of chromosome before the windowstep is completed
				elif position == length_chromosome-1:

					##Infer percent coverage for window
					feature_1_window  = feature_utils.total_window(start_window, position, chromosome, unique_feature_1)	
					feature_1_coverage = float(feature_1_window/count_bases_window)*100

					#Add depth, snp percent of start position
					self.chromosomes_windows[chromosome][start_window] = (feature_1_coverage)
					#Add depth , snp percent  of end position
					self.chromosomes_windows[chromosome][position] = (feature_1_coverage)

				position+=1

###############################################################################################################################################################################
class Plot():

	def __init__(self, fig_name):

		self.fig_name=fig_name

	def plot_features_chromosomes(self, features_chroms, max_len, label_1, species, no_fill):

		fig, axs = plt.subplots(nrows=len(features_chroms), ncols=1, figsize=(14,10), sharey=True, sharex=True, tight_layout=True)
		title = fig.suptitle(label_1+' content along chromosomes'+" - "+species, fontsize=10)
		x_title = fig.supxlabel("Position on chromosome (Mb)")
		y_title = fig.supylabel("Percentage of bases covered in window (%)")
		no=0

		if len(features_chroms.keys())>1:

			for i in features_chroms.keys():
				#Unpack dictionary -> creates lists of values for plotting
				position = []
				y_feat1_coverage = []

				for p, value in sorted(features_chroms[i].items()):
					position.append(p)
					y_feat1_coverage.append(value)

				l1 = axs[no].plot(position, y_feat1_coverage, color="darkorchid", label=label_1)
				if not no_fill:
					axs[no].fill_between(x=position, y1=y_feat1_coverage, color='darkorchid', alpha=.1)
				axs[no].locator_params(axis='y', nbins=4)
				axs[no].set_title(i, fontsize=8)
				axs[no].grid()
				no+=1

		else:
			position = []
			y_feat1_coverage = []

			for p, value in features_chroms[list(features_chroms.keys())[0]].items():
				position.append(p)
				y_feat1_coverage.append(value)

			l1 = axs.plot(position, y_feat1_coverage, color="darkorchid", label=label_1)
			if not no_fill:
				axs.fill_between(x=position, y1=y_feat1_coverage, color='darkorchid', alpha=.1)
			axs.locator_params(axis='y', nbins=4)
			
			axs.set_title(list(features_chroms.keys())[0], fontsize=8)
			axs.grid()

		#adjust intervals of grid depending on the genome size (6 intervals)
		interval=max_len/8
		plt.gca().xaxis.set_major_locator(plt.MultipleLocator(interval))
		current_values = plt.gca().get_xticks().tolist()
		plt.gca().set_xticks(current_values)
		#values in x axis in Mb scale
		plt.gca().set_xticklabels(['{:.1f}'.format(int(x)/1000000) for x in current_values])
		#maximum value on x axis -> length of longest chromosome + 100kb
		plt.xlim(0,max_len+500000)

		if len(features_chroms.keys())>1:
			handles, labels = axs[no-1].get_legend_handles_labels()
		else:
			handles, labels = axs.get_legend_handles_labels()

		lgd = fig.legend(handles, labels, fontsize=8, loc="right", bbox_to_anchor=(1.05, 0.5))
		plt.subplots_adjust(right=1.05)
		fig.savefig(self.fig_name, bbox_extra_artists=(lgd, title, x_title, y_title), bbox_inches='tight')
		plt.show()

####################################################################################################################################################
def main():

	print("\n############################################################################################################################################################")
	print("Estimating feature content (one feature) for chromosome windows and plotting content along the chromosomes of the reference genome...")
	print("############################################################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot content of one feature (e.g., genes or repeats) along the chromosomes using a BED files of genome annotations")
	parser.add_argument("in_assembly", help="Genome assembly fasta used to make the annotations")
	parser.add_argument("feature_1", help="BED file with annotated feature 1 on the reference genome")
	parser.add_argument("feature_1_name", help="Name of feature 1 for plot label without spaces (e.g., coding_sequences)")
	parser.add_argument("out_plot", help="Output plot with coverage, SNPs and repeat content along chromosomes.")
	parser.add_argument("bin_size", help="Size of bins (bp) for plot (average coverage and SNPs percentage plotted for each bin)")
	parser.add_argument("species", help="Species name/sequencing library for title of the plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\" )")
	parser.add_argument('--no_fill', action='store_true', help="Use this flag to indicate no color filling between the lineplot and the x axis")

	parser.usage = 'python3 chromgenomeplot.py feature [positional arguments]'

	#Parsing args
	args = parser.parse_args()
	selected_chrom_to_use=chrom_selection.chromosome_selection(args.sel_chrom)

	#Parsing assembly fasta
	print("Processing Assembly fasta file...\n")
	my_assembly = Assembly_FASTA(filename=args.in_assembly)
	print("Reading Assembly fasta file...\n")
	my_assembly.read_fasta(selected_chrom_to_use)
	print("Finished reading Assembly fasta file...\n"+200*"-"+"\n")

	#Reading feature 1 annotations
	print("Processing BED file with feature annotations...\n")
	my_feature_1 = BEDfile(filename=args.feature_1)
	my_feature_1.read_bed(selected_chrom_to_use)
	print("Finished processing BED file!\n"+200*"-"+"\n")

	#Calculates feature content (percent bases covered) per window along the chromosomes
	print("Extracting average percent feature content for windows of "+str(args.bin_size)+" bp ...")
	my_assembly.calculate_average_bed_values_chromosomes(bin_size=int(args.bin_size), feature_1=my_feature_1.chrom_features)
	print("Finished processing all files!\n"+200*"-"+"\n")
	
	#Return values for plotting (dictionary of dictionaries of tuples)
	#Return dictionary with chromosomes as keys and -> [position -> (feature 1 percent, feature 2 percent)] as value
	#Make plot SNPS and depth along chromosomes
	print("Preparing the plot...")
	my_plot = Plot(fig_name=args.out_plot)
	my_plot.plot_features_chromosomes(features_chroms=my_assembly.chromosomes_windows, max_len=feature_utils.get_max_length(my_assembly.scaffolds_seqs), label_1=args.feature_1_name, species=args.species, no_fill=args.no_fill)
	print("All done!")

########################################################################################################################################################################
if __name__ == '__main__':
	main()

	


