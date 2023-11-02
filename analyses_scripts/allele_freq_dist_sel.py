#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util import chrom_selection, snp_utilities, feature_utils
from allele_freq_dist import VCF_AFD
from allele_freq_chrom import Assembly_FASTA_ALLELE_FREQ
import seaborn as sns
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import argparse

# Checks for potential local variation in dominant allele frequencies (as signal for local aneuploidy)
# Plots 2 distributions on the same figure: 1) Dominant allele frequency distribution when considering entire genome and
# 2) when considering only a specified region a specified region. 

##############################################################################################################################################################
def allelic_freqs_chrom_sel(positions_allelic_freqs, chrom, start, end, bin_size):

	frequencies_chrom_sel = []
	frequencies            = []

	for i in positions_allelic_freqs.keys():

		for p,f in sorted(positions_allelic_freqs[i].items()):
			frequencies.append(f)

			if (p < int(end_size)) or (p > (chrom_lens[i]-int(end_size))):
				frequencies_chrom_ends.append(f)

	return (frequencies, frequencies_chrom_ends)

###############################################################################################################################################################################
class VCF_AFD_SEL(VCF_AFD):

	def __init__(self, filename):

		self.filename = filename
		self.no_snps_uncalled_vcf = 0
		self.no_snps_vcf = 0
		self.no_snps_vcf_filtered = 0
		self.genotype_info_chrom_positions= defaultdict(lambda: [])
		self.positions_uncalled_for_plotting= defaultdict(lambda: [])
		self.positions_snps_filtered= []
		self.homozygous_snps = 0
		self.heterozygous_snps = 0
		self.heterozygous_snps_selection = []
		self.heterozygous_snps_selection_no = 0

	def get_allele_frequencies_selection(self, selected_chrom, start, end, no_snp_filter):

		#Filter SNPS based on separators, biallelic sites, total depth and allelic depth ratio and extract uncalled sites for SNP percent estimation

		for snp in self.genotype_info_chrom_positions.keys():

			#Filter sites based on genotype separator "|" or "\" and return number of alleles for site otherwise exit if unexpected separator
			alleles_for_site=snp_utilities.filter_based_on_separator(self.genotype_info_chrom_positions[snp])
			
			#Check that the calls are bi-allelic
			if snp_utilities.check_biallelic(alleles_for_site):
				
				#Check for uncalled genotypes
				if self.genotype_info_chrom_positions[snp][2][0]=="./." or self.genotype_info_chrom_positions[snp][2][0]==".|.":
					print(f"Warning: Position {self.genotype_info_chrom_positions[snp][1]} on chromosome {self.genotype_info_chrom_positions[snp][0]} does not have called genotypes: GT={self.genotype_info_chrom_positions[snp][2][0]}")

				#check for homozygous (non-variant) sites (could be that the same position has a SNP for another sample in the same VCF so these sites are printed in the combined VCF of SNPs or that there is one SNP that is homozygous but different from the reference)
				elif alleles_for_site[0]==alleles_for_site[1]:
					self.homozygous_snps+=1
					self.no_snps_vcf+=1

				else:
					self.no_snps_vcf+=1
					self.heterozygous_snps+=1
					
					if "/" in self.genotype_info_chrom_positions[snp][2][0]:
						allele_1_number = int(self.genotype_info_chrom_positions[snp][2][0].split("/")[0])
						allele_2_number = int(self.genotype_info_chrom_positions[snp][2][0].split("/")[1])

					elif "|" in self.genotype_info_chrom_positions[snp][2][0]:
						allele_1_number = int(self.genotype_info_chrom_positions[snp][2][0].split("|")[0])
						allele_2_number = int(self.genotype_info_chrom_positions[snp][2][0].split("|")[1])

					ad1= int(self.genotype_info_chrom_positions[snp][2][1].split(",")[allele_1_number])
					ad2= int(self.genotype_info_chrom_positions[snp][2][1].split(",")[allele_2_number])
					
					#Filter based on depth statistics (disable if --no_snp_filter flag is turned on)
					if not no_snp_filter:
						if snp_utilities.check_depth_criteria_het(ad1, ad2):
							self.positions_snps_filtered.append(((max(ad1, ad2))/(ad1+ad2)))
							self.no_snps_vcf_filtered+=1

							#extract AF for selected region
							if self.genotype_info_chrom_positions[snp][0] == selected_chrom and int(self.genotype_info_chrom_positions[snp][1])>=start and int(self.genotype_info_chrom_positions[snp][1])<=end:
								self.heterozygous_snps_selection.append(((max(ad1, ad2))/(ad1+ad2)))
								self.heterozygous_snps_selection_no+=1

							if (max(ad1, ad2)/(ad1+ad2))<0.5:
								print("Error! Some heterozygous SNPs have a frequency of the dominant allele that is smaller than 0.5. Check the format of your VCF file!")
								print(self.genotype_info_chrom_positions[snp])
								exit(1)

					else:
						self.positions_snps_filtered.append(((max(ad1, ad2))/(ad1+ad2)))
						self.no_snps_vcf_filtered+=1

						#extract AF for selected region
						if self.genotype_info_chrom_positions[snp][0] == selected_chrom and int(self.genotype_info_chrom_positions[snp][1])>=start and int(self.genotype_info_chrom_positions[snp][1])<=end:
							self.heterozygous_snps_selection.append(((max(ad1, ad2))/(ad1+ad2)))
							self.heterozygous_snps_selection_no+=1

						if ((max(ad1, ad2))/(ad1+ad2))<0.5:
							print("Error! Some heterozygous SNPs have a frequency of the dominant allele that is smaller than 0.5. Check the format of your VCF file!")
							print(self.genotype_info_chrom_positions[snp])
							exit(1)

			else:
				#Likely multiallelic sites (should not exist for heterozygosity estimates) -> Double check data come if such positions exist in VCF
				print(f"Position {self.genotype_info_chrom_positions[snp][1]} on chromosome {self.genotype_info_chrom_positions[snp][0]} does not have expected genotype format: GT={self.genotype_info_chrom_positions[snp][2][0]}!")
				print("The site likely not biallelic! Please make sure that only bi-allelic sites exist for this sample for accurate SNP-based heterozygosity estimation...\n")
				exit(1)
		
###############################################################################################################################
class Plot():

	def __init__(self, fig_name):

		self.fig_name=fig_name

	def plot_allelic_depth_distribution(self, out_plot, freqs_all, freqs_chrom_sel, bins, bin_size, chromosome, species):

		plt.rcParams["font.family"]= "Arial"

		df_freqs_all        = pd.DataFrame(freqs_all, columns = ["Frequency of dominant allele - full chromosomes"])
		df_freqs_chrom_sel = pd.DataFrame(freqs_chrom_sel, columns = ["Frequency of dominant allele - selected chromosome region"])

		fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 8))
		title = fig.suptitle('Dominant allele frequency distribution - '+species, fontsize=10)
		y_title = fig.supylabel('No. of windows of '+str(bin_size)+"bp", fontsize=10)

		sns.set_style('dark')
		a = sns.histplot(data=df_freqs_all, x="Frequency of dominant allele - full chromosomes", bins=bins, kde = False, color='steelblue', ax=axs[0])
		a.grid(True)
		b = sns.histplot(data=df_freqs_chrom_sel, x="Frequency of dominant allele - selected chromosome region", bins=bins, kde = False, color='olive', ax=axs[1])
		ax=axs[1].set_xlabel( "Frequency of dominant allele - selected chromosome region ("+chromosome+")" , size = 12 )
		b.grid(True)
		current_values = plt.gca().get_yticks().tolist()
		plt.gca().yaxis.set_major_locator(mticker.FixedLocator(current_values))
		plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
		plt.tight_layout()
		plt.savefig(out_plot)
		plt.show()

####################################################################################################################################################
def main():

	print("\n######################################################################################################")
	print("Plotting histogram of dominant allele frequencies from a set of SNPs for a selected chromosome region...")
	print("######################################################################################################\n")


	parser=argparse.ArgumentParser(description="Plot dominant allele frequency distribution for 1) the whome genome and 2) only for selected region")
	parser.add_argument("in_vcf", help="VCF output file of GATK that has been filtered to contain only SNPs. It may contain calls for multiple samples (each on different column).")
	parser.add_argument("column_sample_in_vcf", help="Column no. of genotype info for sample in VCF (starting from 0). Assumes the following format: GT:AD:DP (e.g. 0/1:12,6:18)")
	parser.add_argument("out_plot", help="Output plot name with dominant allele frequency distributions.")
	parser.add_argument("bin_size", help="Bin size (bp) of windows that the average frequency of the dominant allele is inferred")
	parser.add_argument("no_bins", help="No. bins for distribution (histogram) plot")
	parser.add_argument("start_chr", help="Start position of chromosome region for which the frequency of dominant alleles should be inferred")
	parser.add_argument("end_chr", help="End position of chromosome region for which the frequency of dominant alleles should be inferred")
	parser.add_argument("chromosome", help="Selected chromosome")
	parser.add_argument("species", help="Species name/sequencing library for title of the plot")
	parser.add_argument('--no_snp_filter', action='store_true', help="Use this flag if you don't want to filter the SNPs based on total allelic depth for position and allelic depth ratio (heterozygous SNPs)")

	parser.usage='python3 chromgenomeplot.py allele_freq_dist_sel [positional arguments]' 

	#Parsing args
	args=parser.parse_args()
	
	my_vcf = VCF_AFD_SEL(filename=args.in_vcf)
	#Read VCF file and calculate snps for selected sample (column)
	my_vcf.read_extract_snp_variants_for_sample(args.column_sample_in_vcf, selected_chrom_to_use=0)
	print("Total no. of SNPs and uncalled sites: "+str(my_vcf.no_snps_uncalled_vcf)+"\n")
	print("Filtering SNPs from VCF file...\n")
	my_vcf.get_allele_frequencies_selection(selected_chrom=args.chromosome, start=int(args.start_chr), end=int(args.end_chr), no_snp_filter=args.no_snp_filter)
	if args.no_snp_filter:
		print("Filtering of SNPs is disabled...\n")
	print("Total no. of SNPs: "+str(my_vcf.no_snps_vcf))
	print("Total no. of homozygous SNPs: "+str(my_vcf.homozygous_snps))
	print("Total no. of heterozygous SNPs: "+str(my_vcf.heterozygous_snps)+"\n")
	
	if not args.no_snp_filter:
		#Extracting filter snps and uncalled positions
		print("\nTotal no. of heterozygous SNPs that pass the filtering criteria (a) >=20 total depth for position and (2) 0.20 <= allelic depth ratios <= 0.80: "+str(my_vcf.no_snps_vcf_filtered)+"\n")
		print("Total no. of heterozygous SNPs for selected region: "+str(my_vcf.heterozygous_snps_selection_no)+"\n")
		print("Percent of heterozygous SNPs that do not pass the filtering criteria: {percent_removed: .2f}%\n".format(percent_removed=((my_vcf.heterozygous_snps-my_vcf.no_snps_vcf_filtered)/my_vcf.heterozygous_snps)*100))
	else:
		print("No filtering of heterozygous SNPs performed: {percent_removed: .2f}%  of heterozygous SNPs will be used (n={no_used})".format(percent_removed=((my_vcf.no_snps_vcf_filtered)/my_vcf.heterozygous_snps)*100, no_used=my_vcf.no_snps_vcf_filtered))
		print("Total no. of heterozygous SNPs for selected region: "+str(my_vcf.heterozygous_snps_selection_no)+"\n")
	print("Finished processing VCF file!\n"+200*"-"+"\n")

	print("Extracting allele frequencies for target chromosomal region ...\n")
	

	#Make plot SNPS and depth along chromosomes
	print("Preparing the plot...")
	my_plot = Plot(fig_name=args.out_plot)
	my_plot.plot_allelic_depth_distribution(out_plot=args.out_plot, freqs_all=my_vcf.positions_snps_filtered, freqs_chrom_sel=my_vcf.heterozygous_snps_selection, bins=int(args.no_bins), chromosome=args.chromosome, bin_size=int(args.bin_size), species=args.species)
	print("All done!")

########################################################################################################################################################################
if __name__ == '__main__':
	main()

	


