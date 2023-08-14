#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from cov_snps_number import VCF
from feature import Plot
from snps_number import Assembly_FASTA_SNP_NUMBER
from util import chrom_selection, snp_utilities, feature_utils
from collections import defaultdict
import matplotlib.pyplot as plt
import argparse
import os

# Plots snp number for windows along the selected chromosomes
###############################################################################################################################
class Assembly_FASTA_SNP_PERCENT(Assembly_FASTA_SNP_NUMBER):

	def calculate_average_snp_values_chromosomes(self, bin_size, snps):
		'''
		Provide as input window size, snps
		'''
		for chromosome in self.scaffolds_seqs.keys():

			print("Working on chromosome: "+chromosome+" ...")
			#Coverage-specific values
			start_window = 0
			count_bases_window = 0
			length_chromosome = len(self.scaffolds_seqs[chromosome])
			position=0

			for nuc in self.scaffolds_seqs[chromosome]:
				count_bases_window += 1

				if (position+1)%bin_size==0:

					#Infer percent coverage for window
					snps_window = snp_utilities.snps_total_window(start_window, position, chromosome, snps)
					percentage=(snps_window/bin_size)*100

					#Add depth, snp percent of start position
					self.chromosomes_windows[chromosome][start_window] = (percentage)
					#Add depth , snp percent  of end position
					self.chromosomes_windows[chromosome][position] = (percentage)
					start_window= position+1
					count_bases_window = 0

				#if reaching the end of chromosome before the windowstep is completed
				elif position == length_chromosome-1:

					##Infer percent coverage for window
					snps_window  = snp_utilities.snps_total_window(start_window, position, chromosome, snps)
					percentage=(snps_window/bin_size)*100	

					#Add depth, snp percent of start position
					self.chromosomes_windows[chromosome][start_window] = (percentage)
					#Add depth , snp percent  of end position
					self.chromosomes_windows[chromosome][position] = (percentage)

				position+=1

####################################################################################################################################################
class Plot_SNP_NUMBER(Plot):

	def plot_snps_chromosomes(self, snps_chroms, max_len, species, no_fill):

		fig, axs = plt.subplots(nrows=len(snps_chroms), ncols=1, figsize=(14,10), sharey=True, sharex=True, tight_layout=True)
		title = fig.suptitle('SNP percentage in chromosome windows'+" - "+species, fontsize=10)
		x_title = fig.supxlabel("Position on chromosome (Mb)")
		y_title = fig.supylabel("Percentage of SNPs in chromosome windows (%)")
		no=0

		if len(snps_chroms.keys())>1:

			for i in snps_chroms.keys():
				#Unpack dictionary -> creates lists of values for plotting
				position = []
				y_feat1_coverage = []

				for p, value in sorted(snps_chroms[i].items()):
					position.append(p)
					y_feat1_coverage.append(value)

				l1 = axs[no].plot(position, y_feat1_coverage, color="firebrick", label="percentage of SNPs")
				if not no_fill:
					axs[no].fill_between(x=position, y1=y_feat1_coverage, color='red', alpha=.1)
				axs[no].locator_params(axis='y', nbins=4)
				axs[no].set_title(i, fontsize=8)
				axs[no].grid()
				no+=1

		else:
			position = []
			y_feat1_coverage = []

			for p, value in snps_chroms[list(snps_chroms.keys())[0]].items():
				position.append(p)
				y_feat1_coverage.append(value)

			l1 = axs.plot(position, y_feat1_coverage, color="firebrick", label="percentage of SNPs")
			if not no_fill:
				axs.fill_between(x=position, y1=y_feat1_coverage, color='red', alpha=.1)
			axs.locator_params(axis='y', nbins=4)
			
			axs.set_title(list(snps_chroms.keys())[0], fontsize=8)
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

		if len(snps_chroms.keys())>1:
			handles, labels = axs[no-1].get_legend_handles_labels()
		else:
			handles, labels = axs.get_legend_handles_labels()

		lgd = fig.legend(handles, labels, fontsize=8, loc="right", bbox_to_anchor=(1.05, 0.5))
		plt.subplots_adjust(right=1.05)
		fig.savefig(self.fig_name, bbox_extra_artists=(lgd, title, x_title, y_title), bbox_inches='tight')
		plt.show()


def main():

	print("\n###############################################################################################################################################")
	print("Estimating percent of SNPs ( % bases of homozygous + heterozygous SNPs) for chromosome windows and plotting them along the chromosomes of the reference genome ...")
	print("##################################################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot SNP percentages for chromosome windows and plotting along chromosomes of the reference genome")
	parser.add_argument("in_assembly", help="Genome assembly fasta")
	parser.add_argument("in_vcf", help="VCF output file of GATK4 that has been filtered to contain only SNPs. It may contain calls for multiple samples (each on different column).")
	parser.add_argument("column_sample_in_vcf", help="Column no. of genotype information for specified sample in VCF (starting from 0). Assumes the following format: GT:AD:DP (e.g. 0/1:12,6:18)")
	parser.add_argument("out_plot", help="Output plot name with coverage and SNPs plotted along chromosomes")
	parser.add_argument("bin_size", help="Size of bins in bp (SNPs' percentage plotted for each bin)")
	parser.add_argument("species", help="Species name/sequencing library for title of the plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\")")
	parser.add_argument('--no_fill', action='store_true', help="Use this flag to indicate no color filling between the lineplot and the x axis")
	parser.add_argument('--no_snp_filter', action='store_true', help="Use this flag if you don't want to filter the SNPs based on total allelic depth for position")

	parser.usage = 'python3 chromgenomeplot.py snps_percent [positional arguments]'

	#Parsing args
	args = parser.parse_args()
	selected_chrom_to_use=chrom_selection.chromosome_selection(args.sel_chrom)

	#Processing VCF
	print("Processing VCF file...\n")
	#Assign VCF file as instance of VCF class
	my_vcf = VCF(filename=args.in_vcf)
	#Read VCF file and calculate snps for selected sample (column)
	my_vcf.read_extract_snp_variants_for_sample(args.column_sample_in_vcf, selected_chrom_to_use)
	print("Total no. of SNPs and uncalled sites: "+str(my_vcf.no_snps_uncalled_vcf)+"\n")
	print("Filtering SNPs from VCF file...\n")
	if args.no_snp_filter:
		print("Filtering of SNPs is disabled...\n")
	my_vcf.filter_snps(args.no_snp_filter)
	print("Total no. of SNPs: "+str(my_vcf.no_snps_vcf))
	
	if not args.no_snp_filter:
		#Extracting filter snps and uncalled positions
		print("\nTotal no. of SNPs that pass the filtering criteria (a) >=20 total depth for position: "+str(my_vcf.no_snps_vcf_filtered)+"\n")
		print("Percent of SNPs that do not pass the filtering criteria: {percent_removed: .2f}%\n".format(percent_removed=((my_vcf.no_snps_vcf-my_vcf.no_snps_vcf_filtered)/my_vcf.no_snps_vcf)*100))
	else:
		print("No filtering of SNPs performed: {percent_removed: .2f}%  of SNPs will be used (n={no_used})".format(percent_removed=((my_vcf.no_snps_vcf_filtered)/my_vcf.no_snps_vcf)*100, no_used=my_vcf.no_snps_vcf_filtered))
	print("Finished processing VCF file!\n"+200*"-"+"\n")

	#Parsing assembly fasta
	print("Processing Assembly fasta file...\n")
	my_assembly = Assembly_FASTA_SNP_PERCENT(filename=args.in_assembly)
	print("Reading Assembly fasta file...\n")
	my_assembly.read_fasta(selected_chrom_to_use)
	print("Finished reading Assembly fasta file...\n"+200*"-"+"\n")

	#Calculates total SNP number per window along the chromosomes
	print("Extracting total SNP number for windows of "+str(args.bin_size)+"bp...")
	my_assembly.calculate_average_snp_values_chromosomes(bin_size=int(args.bin_size), snps=my_vcf.positions_snps_filtered)
	print("Finished processing all files!\n"+200*"-"+"\n")
		
	#Make plot SNPS numbers along chromosomes
	print("Preparing the plot...")
	my_plot = Plot_SNP_NUMBER(fig_name=args.out_plot)
	my_plot.plot_snps_chromosomes(snps_chroms=my_assembly.chromosomes_windows, max_len=feature_utils.get_max_length(my_assembly.scaffolds_seqs), species=args.species, no_fill=args.no_fill)
	print("All done!\n"+200*"-"+"\n")

###############################################
if __name__ == '__main__':
	main()