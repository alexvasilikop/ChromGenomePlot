#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util import chrom_selection, snp_utilities, feature_utils
from cov_het import VCF
from feature import Assembly_FASTA
from snps_number import Plot_SNP_NUMBER
from collections import defaultdict
import matplotlib.pyplot as plt
import argparse

#Plots dominant allele frequency distribution along chromosomes of the reference genome
#############################################################################################""
def major_allelic_depth_ratio_average(start_window, p, chromosome, snps):

	#count snps for window
	no_snps_window=0
	total_freq=0

	for (i, f) in snps[chromosome]:
		if (i>=start_window) and (i<=p):
			no_snps_window+=1
			total_freq+=f

	if no_snps_window> 0:
		#print("Average freq of dominant allele"
		return float((total_freq/no_snps_window))

	else:
		print("Warning! Some windows contain no heterozygous SNPs. It is recommended to increase window size...")
		print("A default value of 1.0 (frequency of dominant allele) will be printed for windows without SNPs")
		return None

###############################################################################################################################
class Assembly_FASTA_ALLELE_FREQ(Assembly_FASTA):

	def __init__(self, filename):
		super().__init__(filename)
		self.chromosomes_windows_allelic_freq = defaultdict(lambda: {})

	def calculate_average_allelic_frequencies(self, bin_size, snps):

		for chromosome in self.scaffolds_seqs.keys():

			print("Working on chromosome: "+chromosome+" ...")
			#Coverage-specific values
			start_window = 1
			count_bases_window = 0
			length_chromosome = len(self.scaffolds_seqs[chromosome])
			p=0

			for nuc in self.scaffolds_seqs[chromosome]:

				count_bases_window += 1

				if (p+1)%bin_size==0:

					average_allelic_freq = major_allelic_depth_ratio_average(start_window, p+1, chromosome, snps)
					#add snp percent of start position
					self.chromosomes_windows_allelic_freq[chromosome][start_window] = (average_allelic_freq)
					#add snp percent of end position
					self.chromosomes_windows_allelic_freq[chromosome][p+1] = (average_allelic_freq)
					start_window= p+2
					count_bases_window = 0

				#if reaching the end of chromosome before the windowstep is completed
				elif p == length_chromosome-1:

					average_allelic_freq = major_allelic_depth_ratio_average(start_window, p+1, chromosome, snps)

					self.chromosomes_windows_allelic_freq[chromosome][start_window] = (average_allelic_freq)
					self.chromosomes_windows_allelic_freq[chromosome][p+1] = (average_allelic_freq)

				p+=1

#####################################################################################################################################################################
class VCF_AF(VCF):

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

					#store positions with uncalled genotypes for each chromosome to exclude them from the average heterozygosity estimation for each window later
					self.positions_uncalled_for_plotting[self.genotype_info_chrom_positions[snp][0]].append(int(self.genotype_info_chrom_positions[snp][1]))

				#check for homozygous (non-variant) sites (could be that the same position has a SNP for another sample in the same VCF so these sites are printed in the combined VCF of SNPs or that there is one SNP that is homozygous but different from the reference)
				elif alleles_for_site[0]==alleles_for_site[1]:
					self.homozygous_snps+=1
					self.no_snps_vcf+=1

				else:
					self.no_snps_vcf+=1
					self.heterozygous_snps+=1
					self.no_snps_uncalled_vcf+=1
					
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
							self.positions_snps_filtered[self.genotype_info_chrom_positions[snp][0]].append((int(self.genotype_info_chrom_positions[snp][1]),((max(ad1, ad2))/(ad1+ad2))))
							self.no_het_snps_vcf_filtered+=1

					else:
						self.positions_snps_filtered[self.genotype_info_chrom_positions[snp][0]].append((int(self.genotype_info_chrom_positions[snp][1]),((max(ad1, ad2))/(ad1+ad2))))
						self.no_het_snps_vcf_filtered+=1

			else:
				#Likely multiallelic sites (should not exist for heterozygosity estimates) -> Double check data come if such positions exist in VCF
				print(f"Position {self.genotype_info_chrom_positions[snp][1]} on chromosome {self.genotype_info_chrom_positions[snp][0]} does not have expected genotype format: GT={self.genotype_info_chrom_positions[snp][2][0]}!")
				print("The site likely not biallelic! Please make sure that only bi-allelic sites exist for this sample for accurate SNP-based heterozygosity estimation...\n")
				exit(1)

###############################################################################################################################################################################
class Plot_ALLELE_FREQUENCY(Plot_SNP_NUMBER):

	def plot_allele_freq_chromosomes(self, allele_freq, max_len, species, no_fill):

		plt.rcParams["font.family"]= "Arial"

		fig, axs = plt.subplots(nrows=len(allele_freq), ncols=1, figsize=(14,10), sharey=True, sharex=True, tight_layout=True)
		title = fig.suptitle('Dominant allele frequency along chromosomes'+" - "+species, fontsize=10)
		x_title = fig.supxlabel("Position on chromosome (Mb)")
		y_title = fig.supylabel("Dominant allele frequency along chromosomes")
		no=0

		if len(allele_freq.keys())>1:

			for i in allele_freq.keys():
				#Unpack dictionary -> creates lists of values for plotting
				position = []
				y_feat1_coverage = []

				for p, value in sorted(allele_freq[i].items()):
					position.append(p)
					y_feat1_coverage.append(value)

				l1 = axs[no].plot(position, y_feat1_coverage, color="olive", label="allele graquency")
				if not no_fill:
					axs[no].fill_between(x=position, y1=y_feat1_coverage, color='olive', alpha=.1)
				axs[no].locator_params(axis='y', nbins=4)
				axs[no].set_title(i, fontsize=8)
				axs[no].grid()
				no+=1

		else:
			position = []
			y_feat1_coverage = []

			for p, value in allele_freq[list(allele_freq.keys())[0]].items():
				position.append(p)
				y_feat1_coverage.append(value)

			l1 = axs.plot(position, y_feat1_coverage, color="olive", label="allele graquency")
			if not no_fill:
				axs.fill_between(x=position, y1=y_feat1_coverage, color='olive', alpha=.1)
			axs.locator_params(axis='y', nbins=4)
			axs.set_title(list(allele_freq.keys())[0], fontsize=8)
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

		if len(allele_freq.keys())>1:
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
	print("Estimating dominant allele frequency for chromosome windows and plotting them along the chromosomes of the reference genome...")
	print("############################################################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot dominant allele frequency along chromosomes from genome assembly fasta file and VCF file")
	parser.add_argument("in_assembly", help="Genome assembly file in fasta format")
	parser.add_argument("in_vcf", help="VCF output file of GATK4 that has been filtered to contain only SNPs. It may contain calls for multiple samples (each on a different column).")
	parser.add_argument("column_sample_in_vcf", help="Column number of genotype info for specified sample in VCF (starting from 0 for first column). Assumes the following format: GT:AD:DP (e.g. 0/1:12,6:18)")
	parser.add_argument("out_plot", help="Output plot with average dominant allele frequency along chromosomes.")
	parser.add_argument("bin_size", help="Size of bins (bp) for plot (average frequency of dominant alleles plotted for each bin)")
	parser.add_argument("species", help="Species name/sequencing library for title of the plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\" )")
	parser.add_argument('--no_fill', action='store_true', help="Use this flag to indicate no color filling between the lineplot and the x axis")
	parser.add_argument('--no_snp_filter', action='store_true', help="Use this flag if you don't want to filter the SNPs based on total allelic depth for position and allelic depth ratio (for heterozygous SNPs)")

	parser.usage = 'python3 chromgenomeplot.py allele_freq_chrom [positional arguments]'

	#Parsing args
	args = parser.parse_args()
	selected_chrom_to_use=chrom_selection.chromosome_selection(args.sel_chrom)

	#Processing VCF
	print("Processing VCF file...\n")
	#Assign VCF file as instance of VCF class
	my_vcf = VCF_AF(filename=args.in_vcf)
	#Read VCF file and calculate snps for selected sample (column)
	my_vcf.read_extract_snp_variants_for_sample(args.column_sample_in_vcf, selected_chrom_to_use)
	print("Total no. of SNPs and uncalled sites: "+str(my_vcf.no_snps_uncalled_vcf)+"\n")
	print("Filtering SNPs from VCF file...\n")
	if args.no_snp_filter:
		print("Filtering of SNPs is disabled...\n")
	my_vcf.filter_snps(args.no_snp_filter)
	print("Total no. of SNPs: "+str(my_vcf.no_snps_vcf))
	print("Total no. of homozygous SNPs: "+str(my_vcf.homozygous_snps))
	print("Total no. of heterozygous SNPs: "+str(my_vcf.heterozygous_snps))
	
	if not args.no_snp_filter:
		#Extracting filter snps and uncalled positions
		print("\nTotal no. of heterozygous SNPs that pass the filtering criteria (a) >=20 total depth for position and (2) 0.20 <= allelic depth ratios <= 0.80: "+str(my_vcf.no_het_snps_vcf_filtered)+"\n")
		print("Percent of heterozygous SNPs that do not pass the filtering criteria: {percent_removed: .2f}%\n".format(percent_removed=((my_vcf.heterozygous_snps-my_vcf.no_het_snps_vcf_filtered)/my_vcf.heterozygous_snps)*100))
	else:
		print("No filtering of heterozygous SNPs performed: {percent_removed: .2f}%  of heterozygous SNPs will be used (n={no_used})".format(percent_removed=((my_vcf.no_het_snps_vcf_filtered)/my_vcf.heterozygous_snps)*100, no_used=my_vcf.no_het_snps_vcf_filtered))
	print("Finished processing VCF file!\n"+200*"-"+"\n")

	#Parsing assembly fasta
	print("Processing Assembly fasta file...\n")
	my_assembly = Assembly_FASTA_ALLELE_FREQ(filename=args.in_assembly)
	print("Reading Assembly fasta file...\n")
	my_assembly.read_fasta(selected_chrom_to_use)
	print("Finished reading Assembly fasta file...\n"+200*"-"+"\n")

	#Calculates total SNP number per window along the chromosomes
	print("Extracting total SNP number for windows of "+str(args.bin_size)+" bp ...")
	my_assembly.calculate_average_allelic_frequencies(bin_size=int(args.bin_size), snps=my_vcf.positions_snps_filtered)
	print("Finished processing all files!\n"+200*"-"+"\n")
		
	#Make plot SNPS numbers along chromosomes
	print("Preparing the plot...")
	my_plot = Plot_ALLELE_FREQUENCY(fig_name=args.out_plot)
	my_plot.plot_allele_freq_chromosomes(allele_freq=my_assembly.chromosomes_windows_allelic_freq, \
					     max_len=max([len(seq) for seq in my_assembly.scaffolds_seqs.values()]), \
					     species=args.species, \
					     no_fill=args.no_fill)

	print("Average Heterozygosity (average %  of heterozygous SNPs in windows):")
	av_het=snp_utilities.average_heterozygosity(snps=my_vcf.positions_snps_filtered, uncalled=my_vcf.positions_uncalled_for_plotting, genome_size=sum(len(v) for v in my_assembly.scaffolds_seqs.values()))
	print(f"{av_het:.2f}%")
	print("All done!\n"+200*"-"+"\n")

########################################################################################################################################################################
if __name__ == '__main__':
	main()
