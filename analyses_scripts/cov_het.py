#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util.normalize_max import normalize
from util import samtools_depths_util, chrom_selection, snp_utilities
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

# Plots coverage depth and SNP-based heterozygosity along chromosomes based on provided 1) BAM (sorted) and 2) VCF file (from GATK4)
# The VCF file (output of GATK GenotypeGVCFS) has to be filtered to include only SNP variants (e.g. with GATK SelectVariants) -> further filtering of SNPs is done by the script
# (total depth min. = 20, allellic depth proportion >=0.20 and <=0.80)
# If depths file has been generated in previous run it is automatically used (samtools depth -> step is skipped)
# Coverage is printed as normalized to the maximum printed heterozygosity value on the y axis
# A maximum value of % heterozygosity needs to be provided for the plot (e.g. a value of 10 means: 10% max. heterozygosity).

###############################################################################################################################
class BAM_alignment():

	def __init__(self, filename):

		self.filename = filename
		# Dataframe (chromosome, position, depth)
		self.chromosomes_and_depths = defaultdict(lambda: [])
		#Dictionary with lengths of chromosomes
		self.chromosome_lengths = defaultdict(lambda: {})
		#Dictionary of dictionaries: chromosome -> (start/end position in window -> (average depth for window, average snps for window)
		self.chromosomes_windows_depths_snps = defaultdict(lambda: {})
		#average coverage depth
		self.average_depth = 0
		self.genome_size = 0

	def read_depths_file(self, depths_file, selected_chrom_to_use):

		depth_total=0
		(depth_total, self.chromosomes_and_depths, self.genome_size)=samtools_depths_util.read_depths_file_util(depths_file, depth_total, selected_chrom_to_use, self.chromosomes_and_depths, self.genome_size)
		self.average_depth=depth_total/self.genome_size

	def read_alignment(self, selected_chrom_to_use):

		#If depths file has been generated in previous run read this file. Otherwise run samtools depth
		depths_file="depths.txt"

		if os.path.isfile(depths_file):
			print("Found existing depths file! Samtools depth function will be skipped...\n")
			print("Using previously inferred depths.txt file...\n")
			self.read_depths_file(depths_file, selected_chrom_to_use)

		else:
			depths_file="depths.txt"
			samtools_depths_util.run_samtools(self.filename, depths_file)
			self.read_depths_file(depths_file, selected_chrom_to_use)

	def calculate_average_coverages_snps_chromosomes(self, bin_size, snps, uncalled):

		depth_total=0

		'''Provide as input window size, snp, and uncalled position dictionaries {chromX: [p1, p2, ..], chromY: [p1, p2, ..]}'''
		#Infer percent heterozygous snps for each window
		for chromosome in self.chromosomes_and_depths:

			print("Working on chromosome: "+chromosome+" ...")
			#Coverage-specific values
			start_window = 1
			depth_total_window = 0
			count_bases_window = 0
			snps_window = 0
			uncalled_window = 0
			self.chromosome_lengths[chromosome] = len(self.chromosomes_and_depths[chromosome])

			for p, d in self.chromosomes_and_depths[chromosome]:

				depth_total_window += d
				depth_total+=d
				count_bases_window += 1

				if p%bin_size==0:

					#infer percent snps
					snps_window = snp_utilities.snps_total_window(start_window, p, chromosome, snps)
					uncalled_window = snp_utilities.uncalled_total_window(start_window, p, chromosome, uncalled)
					average_depth = int(depth_total_window/bin_size)
					percent_snps = float((snps_window/(bin_size-uncalled_window))*100)

					#add depth, snp percent of start position
					self.chromosomes_windows_depths_snps[chromosome][start_window] = (average_depth, percent_snps)
					#add depth , snp percent  of end position
					self.chromosomes_windows_depths_snps[chromosome][p] = (average_depth, percent_snps)
					depth_total_window = 0
					start_window= p+1
					count_bases_window = 0

				#if reaching the end of chromosome before the windowstep is completed
				elif p == self.chromosome_lengths[chromosome]:

					snps_window = snp_utilities.snps_total_window(start_window, p, chromosome, snps)
					uncalled_window = snp_utilities.uncalled_total_window(start_window, p, chromosome, uncalled)
					average_depth = int(depth_total_window/count_bases_window)
					percent_snps = float((snps_window/(count_bases_window-uncalled_window))*100)

					self.chromosomes_windows_depths_snps[chromosome][start_window] = (average_depth, percent_snps)
					self.chromosomes_windows_depths_snps[chromosome][p] = (average_depth, percent_snps)

		#avearge depth of genome
		self.average_depth=int(depth_total/sum(l for l in self.chromosome_lengths.values()))

############################################################################################################################################################################
class VCF():

	def __init__(self, filename):

		self.filename = filename
		self.no_snps_uncalled_vcf = 0
		self.no_snps_vcf = 0
		self.no_het_snps_vcf_filtered = 0
		self.genotype_info_chrom_positions= defaultdict(lambda: [])
		self.positions_uncalled_for_plotting= defaultdict(lambda: [])
		self.positions_snps_filtered= defaultdict(lambda: [])
		self.homozygous_snps= 0
		self.heterozygous_snps = 0

	def read_extract_snp_variants_for_sample(self, column_no, selected_chrom_to_use):
		'''
		read vcf and extract SNPS+uncalled positions and dictionary with positions, genotypes, depths and allelic depths
		'''
		(self.no_snps_uncalled_vcf, self.genotype_info_chrom_positions)=snp_utilities.extract_snp_variants_for_sample_util(self.filename, self.no_snps_uncalled_vcf, self.genotype_info_chrom_positions, column_no, selected_chrom_to_use)

	def filter_snps(self, no_snp_filter, mode="hom"):

		#Filter SNPS based on separators, biallelic sites, total depth and allelic depth ratio and extract uncalled sites for heterozygous SNP percent estimation

		for snp in self.genotype_info_chrom_positions.keys():

			#Filter sites based on genotype separator "|" or "\" and return number of alleles for site otherwise exit if unexpected separator
			alleles_for_site=snp_utilities.filter_based_on_separator(self.genotype_info_chrom_positions[snp])
			
			#Check that the calls are bi-allelic
			if snp_utilities.check_biallelic(alleles_for_site):
				
				#Check for uncalled genotypes
				if self.genotype_info_chrom_positions[snp][2][0]=="./." or self.genotype_info_chrom_positions[snp][2][0]==".|.":
					print(f"Warning: Position {self.genotype_info_chrom_positions[snp][1]} on chromosome {self.genotype_info_chrom_positions[snp][0]} does not have called genotypes: GT={self.genotype_info_chrom_positions[snp][2][0]}")

					#store positions with uncalled genotypes for each chromosome to exclude them from the average heterozygosity estimation for each window later
					self.positions_uncalled_for_plotting[self.genotype_info_chrom_positions[snp][0]].append(int(self.genotype_info_chrom_positions[snp][1]))
					
				#check for homozygous sites (could be that the same position has a SNP for another sample in the same VCF so these sites are printed in the combined VCF of SNPs or that there is one SNP that is homozygous but different from the reference)
				elif alleles_for_site[0]==alleles_for_site[1]:
					self.homozygous_snps+=1
					self.no_snps_vcf+=1

					#if mode is not specified homozygous SNPs are also included
					if mode=="hom":
						allele_1_number=0
						allele_2_number=0

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
							if snp_utilities.check_allelic_depth_criteria(ad1, ad2):
								self.positions_snps_filtered[self.genotype_info_chrom_positions[snp][0]].append(int(self.genotype_info_chrom_positions[snp][1]))
								self.no_het_snps_vcf_filtered+=1

						else:
							self.positions_snps_filtered[self.genotype_info_chrom_positions[snp][0]].append(int(self.genotype_info_chrom_positions[snp][1]))
							self.no_het_snps_vcf_filtered+=1

				else:
					self.no_snps_vcf+=1
					self.heterozygous_snps+=1
					allele_1_number=0
					allele_2_number=0

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
							self.positions_snps_filtered[self.genotype_info_chrom_positions[snp][0]].append(int(self.genotype_info_chrom_positions[snp][1]))
							self.no_het_snps_vcf_filtered+=1

					else:
						self.positions_snps_filtered[self.genotype_info_chrom_positions[snp][0]].append(int(self.genotype_info_chrom_positions[snp][1]))
						self.no_het_snps_vcf_filtered+=1

			else:
				#Likely multiallelic sites (should not exist for heterozygosity estimates) -> Double check data come if such positions exist in VCF
				print(f"Position {self.genotype_info_chrom_positions[snp][1]} on chromosome {self.genotype_info_chrom_positions[snp][0]} does not have expected genotype format: GT={self.genotype_info_chrom_positions[snp][2][0]}!")
				print("The site likely not biallelic! Please make sure that only bi-allelic sites exist for this sample for accurate SNP-based heterozygosity estimation...\n")
				exit(1)

###############################################################################################################################################################################
class Plot():

	def __init__(self, fig_name):

		self.fig_name=fig_name

	def plot_coverage_heterozygosity_chromosomes(self, snps_and_depths, average_depth, max_het, species, chrom_lengths, no_fill, no_avg_cov):

		plt.rcParams["font.family"]= "Arial"

		fig, axs = plt.subplots(nrows=len(snps_and_depths), ncols=1, figsize=(14,10), sharey=True, sharex=True, tight_layout=True)
		title = fig.suptitle('SNP-based heterozygosity and coverage depth along chromosomes'+" - "+species, fontsize=10)
		x_title = fig.supxlabel("Position on chromosome (Mb)")
		y_title = fig.supylabel("Percentage of heterozygous SNPs in chromosome windows (%)")
		no=0
		
		if len(snps_and_depths.keys())>1:
			for i in snps_and_depths.keys():

				#Unpack dictionary -> creates lists of values for plotting
				position = []
				y_depths = []
				y_snps = []

				for p, value in sorted(snps_and_depths[i].items()):
					position.append(p)
					y_depths.append(value[0])
					y_snps.append(value[1])

				#normalize coverage depth values to be plotted based on the percent heterozygosity y axis
				y_depths_norm = normalize(y_depths, max_het, average_depth, no_avg_cov)

				l1 = axs[no].plot(position, y_snps, color="firebrick", label="%  heterozygous SNPs")
				if not no_fill:
					axs[no].fill_between(x=position, y1=y_snps, color='red', alpha=.1)

				l2 = axs[no].plot(position, y_depths_norm, color="steelblue", label="Coverage depth")
				if not no_fill:
					axs[no].fill_between(x=position, y1=y_depths_norm, color='blue', alpha=.1)
				axs[no].locator_params(axis='y', nbins=4)
				axs[no].set_title(i, fontsize=8)
				axs[no].grid()
				no+=1

		else:
			position = []
			y_depths = []
			y_snps = []

			for p, value in sorted(snps_and_depths[list(snps_and_depths.keys())[0]].items()):
				position.append(p)
				y_depths.append(value[0])
				y_snps.append(value[1])

			#normalize coverage depth values to be plotted based on the percent heterozygosity y axis
			y_depths_norm = normalize(y_depths, max_het, average_depth, no_avg_cov)

			l1 = axs.plot(position, y_snps, color="firebrick", label="%  heterozygous SNPs")
			if not no_fill:
				axs.fill_between(x=position, y1=y_snps, color='red', alpha=.1)

			l2 = axs.plot(position, y_depths_norm, color="steelblue", label="Coverage")
			if not no_fill:
				axs.fill_between(x=position, y1=y_depths_norm, color='blue', alpha=.1)
			axs.locator_params(axis='y', nbins=4)
			axs.set_title(list(snps_and_depths.keys())[0], fontsize=8)
			axs.grid()

		#adjust intervals of grid depending on the genome size (6 intervals)
		interval=max(chrom_lengths.values())/8
		plt.gca().xaxis.set_major_locator(plt.MultipleLocator(interval))
		current_values = plt.gca().get_xticks().tolist()
		plt.gca().set_xticks(current_values)
		#values in x axis in Mb scale
		plt.gca().set_xticklabels(['{:.1f}'.format(int(x)/1000000) for x in current_values])
		#maximum value on x axis -> length of longest chromosome + 100kb
		plt.xlim(0,max(chrom_lengths.values())+500000)

		if len(snps_and_depths.keys())>1:
			handles, labels = axs[no-1].get_legend_handles_labels()
		else:
			handles, labels = axs.get_legend_handles_labels()

		lgd = fig.legend(handles, labels, fontsize=8, loc="right", bbox_to_anchor=(1.05, 0.5))
		plt.subplots_adjust(right=1.05)
		fig.savefig(self.fig_name, bbox_extra_artists=(lgd, title, x_title, y_title), bbox_inches='tight')
		plt.show()

####################################################################################################################################################
def main():

	print("\n###############################################################################################################################################")
	print("Estimating average coverage depth and SNP-based heterozygosity for chromosome windows and plotting them along the chromosomes of the reference genome...")
	print("##################################################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot coverage and SNP-based heterozygosity along chromosomes from BAM and VCF files. Requires SAMtools on path")
	parser.add_argument("in_alignment", help="Genome-read alignment file in BAM format (sorted)")
	parser.add_argument("in_vcf", help="VCF output file of GATK4 that has been filtered to contain only SNPs. It may contain calls for multiple samples (each on different column)")
	parser.add_argument("column_sample_in_vcf", help="Column no. of genotype information for specified sample in VCF (starting from 0). Assumes the following format: GT:AD:DP (e.g. 0/1:12,6:18)")
	parser.add_argument("out_plot", help="Output plot name with coverage and SNPs plotted along chromosomes")
	parser.add_argument("max_het", help="Max. heterozygosity (percent heterozygous SNPs) to plot on y axis (for coverage normalization). If the flag --no_avg_cov is used then coverage is not normalized to be comparable across chromosomes but only based on a max_het value")
	parser.add_argument("bin_size", help="Size of bins (bp) for coverage plot (average coverage and SNPs percentage plotted for each bin)")
	parser.add_argument("species", help="Species name/sequencing library for title of the plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\")")
	parser.add_argument('--no_fill', action='store_true', help="Use this flag to indicate no color filling between the lineplot and the x axis")
	parser.add_argument('--no_snp_filter', action='store_true', help="Use this flag if you don't want to filter the SNPs based on total allelic depth for position and allelic depth ratio")
	parser.add_argument('--no_avg_cov', action='store_true', help="Use this flag if you don't want to normalize the coverage by its average so that it is comparable among chromosomes (coverage is only normalized to the max_het value)")

	parser.usage = 'python3 chromgenomeplot.py cov_het [positional arguments]'

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
		print("Filtering of heterozygous SNPs is disabled...\n")
	my_vcf.filter_snps(args.no_snp_filter, mode="het")
	print("Total no. of SNPs: "+str(my_vcf.no_snps_vcf))
	print("Total no. of homozygous SNPs: "+str(my_vcf.homozygous_snps))
	print("Total no. of heterozygous SNPs: "+str(my_vcf.heterozygous_snps)+"\n")
	
	if not args.no_snp_filter:
		#Extracting filter snps and uncalled positions
		print("\nTotal no. of heterozygous SNPs that pass the filtering criteria (a) >=20 total depth for position and (2) 0.20 <= allelic depth ratios <= 0.80: "+str(my_vcf.no_het_snps_vcf_filtered)+"\n")
		print("Percent of heterozygous SNPs that do not pass the filtering criteria: {percent_removed: .2f}%\n".format(percent_removed=((my_vcf.heterozygous_snps-my_vcf.no_het_snps_vcf_filtered)/my_vcf.heterozygous_snps)*100))
	else:
		print("No filtering of heterozygous SNPs performed: {percent_removed: .2f}%  of heterozygous SNPs will be used (n={no_used})".format(percent_removed=((my_vcf.no_het_snps_vcf_filtered)/my_vcf.heterozygous_snps)*100, no_used=my_vcf.no_het_snps_vcf_filtered))
	print("Finished processing VCF file!\n"+200*"-"+"\n")

	#Processing BAM alignment
	#Assign bam file as an instance of the Alignment Class
	print("Processing BAM_alignment file...\n")
	my_alignment=BAM_alignment(filename=args.in_alignment)
	#Read BAM file and calculate per base coverages -> generates depths.txt file
	print("Running samtools depth on BAM alignment file...\n")
	my_alignment.read_alignment(selected_chrom_to_use)
	print("Genome size from depths file (all genome positions): "+str(my_alignment.genome_size)+"\n")

	#Calculates coverages and snps per window along the chromosomes
	print("Extracting average coverages and percent heterozygous SNPs for chromosome windows of "+str(args.bin_size)+"bp...\n")
	my_alignment.calculate_average_coverages_snps_chromosomes(bin_size=int(args.bin_size), snps=my_vcf.positions_snps_filtered, uncalled=my_vcf.positions_uncalled_for_plotting)
	print("Finished processing BAM/depths/snps files!\n"+200*"-"+"\n")
		
	#Return values for plotting (dictionary of dictionaries of tuples)
	#Return dictionary with chromosomes as keys and -> [position -> (coverage, snps)] as value
	print("Returning values for plotting...\n")
	
	#Make plot SNPS and depth along chromosomes
	print("Preparing the plot...\n")
	my_plot = Plot(fig_name=args.out_plot)

	my_plot.plot_coverage_heterozygosity_chromosomes(snps_and_depths=my_alignment.chromosomes_windows_depths_snps, \
		average_depth= my_alignment.average_depth, \
		max_het=int(args.max_het), \
		species=args.species, \
		chrom_lengths=my_alignment.chromosome_lengths, \
		no_fill=args.no_fill, \
		no_avg_cov=args.no_avg_cov)

	print("All done!\n"+200*"-"+"\n")

	#Print average % het SNPs and coverage depth for genome
	print("\n"+"Average coverage depth:")
	print(f"{my_alignment.average_depth:.0f} X")

	print("\n% Heterozygosity of the genome (% of bases with heterozygous SNPs):")
	av_het=snp_utilities.average_heterozygosity(snps=my_vcf.positions_snps_filtered, uncalled=my_vcf.positions_uncalled_for_plotting, genome_size=my_alignment.genome_size)
	print(f"{av_het:.2f}%")

	print("\nAverage percent Heterozygosity (Average % of heterozygous SNPs in chromosome windows):")
	print(f"{np.mean(snp_utilities.average_no_snps_windows(my_alignment.chromosomes_windows_depths_snps)):.2f} %")
	

###############################################
if __name__ == '__main__':
	main()
