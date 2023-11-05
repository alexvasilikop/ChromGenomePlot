#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util import chrom_selection, feature_utils
from collections import defaultdict
import matplotlib.pyplot as plt
import argparse
import re
import numpy as np

###########################################################################################################################################################################################################################
# Plots one feature along chromosomes (e.g. gene content or repetitive content) using CHROM files of annotations (each CHROM annotation must contain sorted annotations by genome position)
##############################################################################################################################################################################################################
class CHROMfile():

	def __init__(self, filename):

		self.filename = filename
		self.chrom_features = defaultdict(lambda: [])

	def read_CHROM(self, selected_chrom_to_use):

		self.chrom_features=feature_utils.CHROM_read_util(self.filename, self.chrom_features, selected_chrom_to_use)

########################################################################################################################################
class Assembly_FASTA():

	def __init__(self, filename):

		self.filename = filename
		self.scaffolds_seqs = defaultdict(lambda: "")
		self.chromosomes_windows_coverage_numbers = defaultdict(lambda: {})

	def read_fasta(self, selected_chrom_to_use):

		self.scaffolds_seqs=feature_utils.fasta_read_util(self.filename, self.scaffolds_seqs, selected_chrom_to_use)

	def calculate_coverage_number_values_chromosomes(self, bin_size, feature_1, feature_1_label, mode_numbers):

		#return percent coverage from windows
		self.chromosomes_windows_coverage_numbers = feature_utils.total_window_numbers_or_coverage(bin_size, self.scaffolds_seqs, feature_1, feature_1_label, mode_numbers)

###############################################################################################################################################################################
class Plot():

	def __init__(self, fig_name):

		self.fig_name=fig_name

	def plot_features_chromosomes(self, features_chroms, max_len, label_1, species, no_fill, bin_size, mode_numbers):

		plt.rcParams["font.family"]= "Arial"

		fig, axs = plt.subplots(nrows=len(features_chroms), ncols=1, figsize=(14,10), sharey=True, sharex=True, tight_layout=True)
		title = fig.suptitle('Content of '+label_1+' along chromosomes/scaffolds'+" - "+species, fontsize=10)
		x_title = fig.supxlabel("Position on chromosome (Mb)")
		if mode_numbers:
			y_title = fig.supylabel("Number of features in windows of size "+str(bin_size)+" bp")
		else:
			y_title = fig.supylabel("Percentage of bases covered in chromosome windows of size "+str(bin_size)+" bp (%)")

		no=0
		if len(features_chroms.keys())>1:
			#for more than one chromosome

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

	parser=argparse.ArgumentParser(description="Plot content of one feature (e.g., genes, coding sequences or repeats) along the chromosomes using a CHROM file of genome annotations")
	parser.add_argument("in_assembly", help="Genome assembly fasta used to make the annotations")
	parser.add_argument("feature_1", help="CHROM coordinate file with annotated feature on the reference genome")
	parser.add_argument("feature_1_name", help="Name of feature for plot label without spaces (e.g., coding_sequences, introns)")
	parser.add_argument("out_plot", help="Output plot with feature content along chromosomes")
	parser.add_argument("bin_size", help="Size of bins (bp) for plot (percent bases covered by feature plotted for each bin)/or number of features if --numbers")
	parser.add_argument("species", help="Species name / sequencing library for title of the plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\" )")
	parser.add_argument('--no_fill', action='store_true', help="Use this flag to indicate no color filling between the lineplot and the x axis")
	parser.add_argument('--numbers', action='store_true', help="Use this flag to indicate numbers of features instead of coverage of bases in chromosome windows")

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
	print("Processing CHROM file with feature annotations...\n")
	my_feature_1 = CHROMfile(filename=args.feature_1)
	#Read CHROM
	my_feature_1.read_CHROM(selected_chrom_to_use)
	print("Finished processing CHROM file!\n"+200*"-"+"\n")

	print("Extracting feature content for windows of "+str(args.bin_size)+" bp ...")
	my_assembly.calculate_coverage_number_values_chromosomes(bin_size=int(args.bin_size), feature_1=my_feature_1.chrom_features, feature_1_label=args.feature_1_name, mode_numbers=args.numbers)
	print("Finished processing all files!\n"+200*"-"+"\n")
	
	#Make plot feature number or coverage along chromosomes
	print("Preparing the plot...")
	my_plot = Plot(fig_name=args.out_plot)
	if not args.numbers:
		my_plot.plot_features_chromosomes(features_chroms=my_assembly.chromosomes_windows_coverage_numbers, \
						  max_len=max([len(seq) for seq in my_assembly.scaffolds_seqs.values()]), \
						  label_1=args.feature_1_name, \
						  species=args.species, \
						  no_fill=args.no_fill, \
						  bin_size=args.bin_size, \
						  mode_numbers=args.numbers)
		print("\n## Average base coverage of feature in each chromosome window: {percent:.2f} % ##\n".format(percent=np.mean(feature_utils.coverage_or_numbers_per_window(my_assembly.chromosomes_windows_coverage_numbers, args.numbers))))

	else:
		my_plot.plot_features_chromosomes(features_chroms=my_assembly.chromosomes_windows_coverage_numbers, \
						  max_len=max([len(seq) for seq in my_assembly.scaffolds_seqs.values()]), \
						  label_1=args.feature_1_name, \
						  species=args.species, \
						  no_fill=args.no_fill,
						  bin_size=args.bin_size, \
						  mode_numbers=args.numbers)
		print("\n## Average number of features in each chromosome window {percent:.2f} ##\n".format(percent=np.mean(feature_utils.coverage_or_numbers_per_window(my_assembly.chromosomes_windows_coverage_numbers, args.numbers))))

	print("All done!")

########################################################################################################################################################################
if __name__ == '__main__':
	main()

	


