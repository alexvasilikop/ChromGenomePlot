#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util import chrom_selection, feature_utils
from feature import CHROMfile, Assembly_FASTA
from feat_2_feat import Plot
from collections import defaultdict
import matplotlib.pyplot as plt
import argparse

###########################################################################################################################################################################################################################
# Plots GC content and feature content along chromosomes using CHROM file of annotations (each bed annotation must contain sorted annotations by genome position annotations)
##############################################################################################################################################################################################################

class Assembly_FASTA_new(Assembly_FASTA):

	def calculate_coverage_number_values_chromosomes(self, bin_size, feature_1, feature_1_label, mode_numbers=False):
		self.chromosomes_windows_coverage_numbers = feature_utils.total_window_numbers_or_coverage(bin_size, self.scaffolds_seqs, feature_1, feature_1_label, mode_numbers)
		#Get intervals with GC content and add them to the feature values (%)
		self.chromosomes_windows_coverage_numbers = feature_utils.get_GC_content_for_windows(self.chromosomes_windows_coverage_numbers, self.scaffolds_seqs, bin_size)

####################################################################################################################################################
def main():

	print("\n############################################################################################################################################################")
	print("Estimating GC and feature content for chromosome windows and plotting along the chromosomes of the reference genome...")
	print("############################################################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot GC and feature content (e.g., genes, repeats) along the chromosomes")
	parser.add_argument("in_assembly", help="Genome assembly fasta used to make the annotations")
	parser.add_argument("feature_1", help="CHROM file with annotated feature 1 on the reference genome")
	parser.add_argument("feature_1_name", help="Name of feature 1 for plot label without spaces (e.g., coding_sequences)")
	parser.add_argument("out_plot", help="Output plot with GC content and feature content along the chromosomes")
	parser.add_argument("bin_size", help="Size of bins (bp) for plot (average coverage and SNPs percentage plotted for each bin)")
	parser.add_argument("species", help="Species name/sequencing library for title of the plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\" )")
	parser.add_argument('--no_fill', action='store_true', help="Use this flag to indicate no color filling between the lineplot and the x axis")

	parser.usage = 'python3 chromgenomeplot.py GC_2_feat [positional arguments]'

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

	#Calculates feature content (percent bases covered) per window along the chromosomes
	print("Extracting average percent "+args.feature_1_name+" and GC content (%) for windows of "+str(args.bin_size)+" bp ...")
	my_assembly.calculate_coverage_number_values_chromosomes(bin_size=int(args.bin_size), feature_1=my_feature_1.chrom_features, feature_1_label=args.feature_1_name)
	print("Finished processing all files!\n"+200*"-"+"\n")
	
	#Plotting
	#Make plot
	print("Preparing the plot...")
	my_plot = Plot(fig_name=args.out_plot)
	my_plot.plot_features_chromosomes(bin_size=args.bin_size, \
		                          features_chroms=my_assembly.chromosomes_windows_coverage_numbers, \
					  max_len=max([len(seq) for seq in my_assembly.scaffolds_seqs.values()]), \
					  label_1=args.feature_1_name, \
					  label_2="GC content (%)", \
				          species=args.species, \
					  no_fill=args.no_fill, \
					  mode_numbers=False)
	print("All done!")

########################################################################################################################################################################
if __name__ == '__main__':
	main()

	


