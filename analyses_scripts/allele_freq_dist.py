#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util import chrom_selection, snp_utilities
from cov_het import VCF
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import argparse

#Plots a histogram (distribution) of dominant allele frequencies from a VCF file.
###############################################################################################################################################################################
class VCF_AFD(VCF):

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

	def filter_snps(self, no_snp_filter):
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
					ad1= int(self.genotype_info_chrom_positions[snp][2][1].split(",")[0])
					ad2= int(self.genotype_info_chrom_positions[snp][2][1].split(",")[1])

					#Filter based on depth statistics (disable if --no_snp_filter flag is turned on)
					if not no_snp_filter:
						if snp_utilities.check_depth_criteria_het(ad1, ad2):
							self.positions_snps_filtered.append(((max(ad1, ad2))/(ad1+ad2)))
							self.no_snps_vcf_filtered+=1

							if (max(ad1, ad2)/(ad1+ad2))<0.5:
								print("Error! Some heterozygous SNPs have a frequency of the dominant allele that is smaller than 0.5. Check the format of your VCF file!")
								print(self.genotype_info_chrom_positions[snp])
								exit(1)

					else:
						self.positions_snps_filtered.append(((max(ad1, ad2))/(ad1+ad2)))
						self.no_snps_vcf_filtered+=1
						if ((max(ad1, ad2))/(ad1+ad2))<0.5:
								print("Error! Some heterozygous SNPs have a frequency of the dominant allele that is smaller than 0.5. Check the format of your VCF file!")
								print(self.genotype_info_chrom_positions[snp])

			else:
				#Likely multiallelic sites (should not exist for heterozygosity estimates) -> Double check data come if such positions exist in VCF
				print(f"Position {self.genotype_info_chrom_positions[snp][1]} on chromosome {self.genotype_info_chrom_positions[snp][0]} does not have expected genotype format: GT={self.genotype_info_chrom_positions[snp][2][0]}!")
				print("The site likely not biallelic! Please make sure that only bi-allelic sites exist for this sample for accurate SNP-based heterozygosity estimation...\n")
				exit(1)

############################################################################################################################################################################"""
class Plot():

	def __init__(self, fig_name):

		self.fig_name=fig_name

	def plot_allelic_depth_dist(self, allelic_depths_frequencies, bins, species, out_plot):

		sns.set_style("darkgrid")
		sns.histplot(data=allelic_depths_frequencies, bins=int(bins), kde = False, color='steelblue')
		x_title = plt.xlabel('Frequency of dominant allele', fontsize=8)
		y_title = plt.ylabel('No. of heterozygous SNPs', fontsize=8)
		title = plt.title('Histogram of dominant allele frequencies - '+species, fontsize=10)
		current_values = plt.gca().get_yticks().tolist()
		plt.gca().yaxis.set_major_locator(mticker.FixedLocator(current_values))
		plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
		plt.grid(True)
		plt.tight_layout()
		plt.savefig(out_plot, bbox_extra_artists=(title, x_title, y_title))
		plt.show()

####################################################################################################################################################
def main():

	print("\n############################################################################################################################################################")
	print("Plotting histogram of dominant allele frequencies from a set of SNPs...")
	print("############################################################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot histogram (distribution plot) of dominant allele frequencies using a set of SNPs from a VCF file")
	parser.add_argument("in_vcf", help="VCF output file of GATK that has been filtered to contain only SNPs. It may contain calls for multiple samples (each on different column).")
	parser.add_argument("column_sample_in_vcf", help="Column number of genotype info for specified sample in VCF (starting from 0 for first column). Assumes the following format: GT:AD:DP (e.g. 0/1:12,6:18)")
	parser.add_argument("out_plot", help="Output plot name (dominant allele frequency distribution plot)")
	parser.add_argument("no_bins", help="No. of bins for allelic depth distribution plot")
	parser.add_argument("species", help="Species name/sequencing library for title of the plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\")")
	parser.add_argument('--no_snp_filter', action='store_true', help="Use this flag if you don't want to filter the SNPs based on total allelic depth for position")

	parser.usage = 'python3 chromgenomeplot.py allele_freq_dist [positional arguments]'
	args=parser.parse_args()
	selected_chrom_to_use=chrom_selection.chromosome_selection(args.sel_chrom)
	
	my_vcf = VCF_AFD(filename=args.in_vcf)
	#Read VCF file and calculate snps for selected sample (column)
	my_vcf.read_extract_snp_variants_for_sample(args.column_sample_in_vcf, selected_chrom_to_use)
	print("Total no. of SNPs and uncalled sites: "+str(my_vcf.no_snps_uncalled_vcf)+"\n")
	print("Filtering SNPs from VCF file...\n")
	if args.no_snp_filter:
		print("Filtering of SNPs is disabled...\n")
	my_vcf.filter_snps(args.no_snp_filter)
	print("Total no. of SNPs: "+str(my_vcf.no_snps_vcf))
	print("Total no. of homozygous SNPs: "+str(my_vcf.homozygous_snps))
	print("Total no. of heterozygous SNPs: "+str(my_vcf.heterozygous_snps)+"\n")
	
	if not args.no_snp_filter:
		#Extracting filter snps and uncalled positions
		print("\nTotal no. of heterozygous SNPs that pass the filtering criteria (a) >=20 total depth for position and (2) 0.20 <= allelic depth ratios <= 0.80: "+str(my_vcf.no_snps_vcf_filtered)+"\n")
		print("Percent of heterozygous SNPs that do not pass the filtering criteria: {percent_removed: .2f}%\n".format(percent_removed=((my_vcf.heterozygous_snps-my_vcf.no_snps_vcf_filtered)/my_vcf.heterozygous_snps)*100))
	else:
		print("No filtering of heterozygous SNPs performed: {percent_removed: .2f}%  of heterozygous SNPs will be used (n={no_used})".format(percent_removed=((my_vcf.no_snps_vcf_filtered)/my_vcf.heterozygous_snps)*100, no_used=my_vcf.no_snps_vcf_filtered))
	print("Finished processing VCF file!\n"+200*"-"+"\n")

	#Plot dominant allele distribution
	print("Preparing the plot...")
	my_plot = Plot(fig_name=args.out_plot)
	my_plot.plot_allelic_depth_dist(allelic_depths_frequencies=my_vcf.positions_snps_filtered, bins= args.no_bins, species=args.species, out_plot=args.out_plot)
	print("All done!")

########################################################################################################################################################################

if __name__ == '__main__':
	main()

	


