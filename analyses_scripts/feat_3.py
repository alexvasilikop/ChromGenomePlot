#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util import chrom_selection, feature_utils
from feature import CHROMfile, Assembly_FASTA, Plot
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

###########################################################################################################################################################################################################################
# Plots two features along chromosomes (e.g. gene content vs repetitive content) using CHROM files of annotations (each CHROM annotation must contain sorted annotations by genome position)
##############################################################################################################################################################################################################

class Assembly_FASTA_new(Assembly_FASTA):

	def __init__(self, filename):
		Assembly_FASTA.__init__(self, filename)
		self.chromosomes_windows_coverage_numbers_feat_1=defaultdict(lambda: {})
		self.chromosomes_windows_coverage_numbers_feat_2=defaultdict(lambda: {})
		self.chromosomes_windows_coverage_numbers_feat_3=defaultdict(lambda: {})

	def calculate_coverage_values_chromosomes(self, bin_size, feature_1, feature_1_label, feature_2, feature_2_label, feature_3, feature_3_label, mode_numbers):
		#return numbers or coverages per window for plotting
		self.chromosomes_windows_coverage_numbers_feat_1 = feature_utils.total_window_numbers_or_coverage(bin_size, self.scaffolds_seqs, feature_1, feature_1_label, mode_numbers)
		self.chromosomes_windows_coverage_numbers_feat_2 = feature_utils.total_window_numbers_or_coverage(bin_size, self.scaffolds_seqs, feature_2, feature_2_label, mode_numbers)
		self.chromosomes_windows_coverage_numbers_feat_3 = feature_utils.total_window_numbers_or_coverage(bin_size, self.scaffolds_seqs, feature_3, feature_3_label, mode_numbers)

		#create data structure for three features
		for c in self.chromosomes_windows_coverage_numbers_feat_2:
			for pos in self.chromosomes_windows_coverage_numbers_feat_2[c]:
				self.chromosomes_windows_coverage_numbers_feat_1[c][pos]=[self.chromosomes_windows_coverage_numbers_feat_1[c][pos], \
				                                                          self.chromosomes_windows_coverage_numbers_feat_2[c][pos], \
				                                                          self.chromosomes_windows_coverage_numbers_feat_3[c][pos]]

###############################################################################################################################################################################
class Plot(Plot):

	def plot_features_chromosomes(self, bin_size, features_chroms, max_len, label_1, label_2, label_3, species, no_fill, mode_numbers):

		plt.rcParams["font.family"]= "Arial"

		fig, axs = plt.subplots(nrows=len(features_chroms), ncols=1, figsize=(14,10), sharey=True, sharex=True, tight_layout=True)
		title = fig.suptitle("Content of features along chromosomes "+label_1+", "+label_2+" and "+label_3+" - "+species, fontsize=10)
		x_title = fig.supxlabel("Position on chromosome (Mb)")

		#y title acccording to mode
		if mode_numbers:
			y_title = fig.supylabel("Number of features in windows of size "+str(bin_size)+" bp")
		else:
			y_title = fig.supylabel("Percentage of bases covered in chromosome windows of size "+str(bin_size)+" bp (%)")

		no=0

		if len(features_chroms.keys())>1:
			for i in features_chroms.keys():

				#Unpack dictionary -> creates lists of values for plotting
				position = []
				y_feat1_coverage = []
				y_feat2_coverage = []
				y_feat3_coverage = []

				for p, value in sorted(features_chroms[i].items()):
					position.append(p)
					y_feat1_coverage.append(value[0])
					y_feat2_coverage.append(value[1])
					y_feat3_coverage.append(value[2])

				l1 = axs[no].plot(position, y_feat1_coverage, color="darkorchid", label=label_1)
				if not no_fill:
					axs[no].fill_between(x=position, y1=y_feat1_coverage, color='darkorchid', alpha=.1)
				axs[no].locator_params(axis='y', nbins=4)
			
				l2 = axs[no].plot(position, y_feat2_coverage, color="darkgoldenrod", label=label_2)
				if not no_fill:
					axs[no].fill_between(x=position, y1=y_feat2_coverage, color='darkgoldenrod', alpha=.1)
				axs[no].locator_params(axis='y', nbins=4)

				l3 = axs[no].plot(position, y_feat3_coverage, color="royalblue", label=label_3)
				if not no_fill:
					axs[no].fill_between(x=position, y1=y_feat3_coverage, color='royalblue', alpha=.1)
				axs[no].locator_params(axis='y', nbins=4)

				axs[no].set_title(i, fontsize=8)
				axs[no].grid()
				no+=1

		else:
			position = []
			y_feat1_coverage = []
			y_feat2_coverage = []
			y_feat3_coverage = []
			
			for p, value in features_chroms[list(features_chroms.keys())[0]].items():
				position.append(p)
				y_feat1_coverage.append(value[0])
				y_feat2_coverage.append(value[1])
				y_feat3_coverage.append(value[2])

			l1 = axs.plot(position, y_feat1_coverage, color="darkorchid", label=label_1)
			if not no_fill:
				axs.fill_between(x=position, y1=y_feat1_coverage, color='darkorchid', alpha=.1)
			axs.locator_params(axis='y', nbins=4)
			
			l2 = axs.plot(position, y_feat2_coverage, color="darkgoldenrod", label=label_2)
			if not no_fill:
				axs.fill_between(x=position, y1=y_feat2_coverage, color='darkgoldenrod', alpha=.1)
			axs.locator_params(axis='y', nbins=4)

			l3 = axs.plot(position, y_feat3_coverage, color="royalblue", label=label_3)
			if not no_fill:
				axs.fill_between(x=position, y1=y_feat3_coverage, color='royalblue', alpha=.1)
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
	print("Estimating feature content (3 features) for chromosome windows and plotting fearure content along the chromosomes of the reference genome...")
	print("############################################################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot content of three features (e.g., genes, repeats and introns) along the chromosomes using CHROM files of genome annotations")
	parser.add_argument("in_assembly", help="Genome assembly fasta used to make the annotations.")
	parser.add_argument("feature_1", help="CHROM file with annotated feature 1 on the reference genome")
	parser.add_argument("feature_2", help="CHROM file with annotated feature 2 on the reference genome")
	parser.add_argument("feature_3", help="CHROM file with annotated feature 3 on the reference genome")
	parser.add_argument("feature_1_name", help="Name of feature 1 for plot label without spaces (e.g., coding_sequences)")
	parser.add_argument("feature_2_name", help="Name of feature 2 for plot label without spaces (e.g., repeat_content)")
	parser.add_argument("feature_3_name", help="Name of feature 3 for plot label without spaces (e.g., introns)")
	parser.add_argument("out_plot", help="Output plot with feature 1 and feature 2 content along chromosomes")
	parser.add_argument("bin_size", help="Size of bins (bp) for plot")
	parser.add_argument("species", help="Species name/sequencing library for title of the plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\" )")
	parser.add_argument('--no_fill', action='store_true', help="Use this flag to indicate no color filling between the lineplot and the x axis")
	parser.add_argument('--numbers', action='store_true', help="Use this flag to indicate numbers of features instead of coverage of bases in chromosome windows")

	parser.usage = 'python3 chromgenomeplot.py feat_3 [positional arguments]'

	#Parsing args
	args = parser.parse_args()
	selected_chromosomes=chrom_selection.chromosome_selection(args.sel_chrom)

	#Parsing assembly fasta
	print("\nProcessing Assembly fasta file...")
	my_assembly = Assembly_FASTA_new(filename=args.in_assembly)
	print("Reading Assembly fasta file...")
	my_assembly.read_fasta(selected_chromosomes)
	print("Finished reading Assembly fasta file...\n"+200*"-"+"\n")

	#Reading feature 1 annotations
	print("Processing CHROM file with feature 1 annotations...")
	my_feature_1 = CHROMfile(filename=args.feature_1)
	my_feature_1.read_CHROM(selected_chromosomes)
	print("Finished processing CHROM file!\n"+200*"-"+"\n")

	#Reading feature 2 annotations
	print("Processing CHROM file with feature 2 annotations...")
	my_feature_2 = CHROMfile(filename=args.feature_2)
	my_feature_2.read_CHROM(selected_chromosomes)
	print("Finished processing CHROM file!\n"+200*"-"+"\n")

	#Reading feature 3 annotations
	print("Processing CHROM file with feature 3 annotations...")
	my_feature_3 = CHROMfile(filename=args.feature_3)
	my_feature_3.read_CHROM(selected_chromosomes)
	print("Finished processing CHROM file!\n"+200*"-"+"\n")

	#Calculates feature content (percent bases covered or numbers) per window along the chromosomes
	print("Extracting "+args.feature_1_name+", "+args.feature_2_name+" and "+args.feature_3_name+" content for windows of "+str(args.bin_size)+" bp ...")
	my_assembly.calculate_coverage_values_chromosomes(bin_size=int(args.bin_size), \
		                                              feature_1=my_feature_1.chrom_features, feature_1_label=args.feature_1_name, \
		                                              feature_2=my_feature_2.chrom_features, feature_2_label=args.feature_2_name, \
													  feature_3=my_feature_3.chrom_features, feature_3_label=args.feature_3_name, \
													  mode_numbers=args.numbers)
	print("Finished processing all files!\n"+200*"-"+"\n")
	
	#Make plot
	my_plot = Plot(fig_name=args.out_plot)
	my_plot.plot_features_chromosomes(bin_size=args.bin_size, \
		                              features_chroms=my_assembly.chromosomes_windows_coverage_numbers_feat_1, \
								      max_len=max([len(seq) for seq in my_assembly.scaffolds_seqs.values()]), \
									  label_1=args.feature_1_name, \
									  label_2=args.feature_2_name, \
									  label_3=args.feature_3_name, \
									  species=args.species, \
									  no_fill=args.no_fill, \
									  mode_numbers=args.numbers)

	print("All done!")

########################################################################################################################################################################
if __name__ == '__main__':
	main()