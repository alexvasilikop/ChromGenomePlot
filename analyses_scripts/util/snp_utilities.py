#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Utilities for processing VCF files and performing reading, filtering and checks of SNPs
#############################################################################################""
def snps_total_window(start_window, p, chromosome, snps):

	#count snps for window
	no_snps_window=0
	for i in snps[chromosome]:
		if (i>=start_window) and (i<=p):
			no_snps_window+=1

	return int(no_snps_window)

###################################################################################
def uncalled_total_window(start_window, p, chromosome, uncalled):

	#count uncalled sites for window
	no_uncalled_window=0
	for i in uncalled[chromosome]:
		if (i>=start_window) and (i<=p):
			no_uncalled_window+=1

	return int(no_uncalled_window)

##############################################################################
def average_heterozygosity(snps, uncalled, genome_size):

	#measures average snp-based heterozygosity not average of windows but average heterozygous snps divided by genomesize-uncalled sites
	'''Input: dictionary of chrom->list_of_snp_positions, dictionary of chrom->list_of_uncalled_positions, genome size'''
	no_snps=0
	no_uncalled=0

	for chrom in sorted(snps.keys()):
		for snp in snps[chrom]:
			no_snps+=1

	for chrom in uncalled.keys():
		for uncalled_site in uncalled[chrom]:
			no_uncalled+=1

	return (no_snps/(genome_size-no_uncalled)) *100

#####################################################################
def average_no_snps_windows(chromosomes_positions_values, cov=True):
	#returns snps numbers / percentages for all windows as a list
	values=[]

	if cov:
		for c in chromosomes_positions_values:
			index=0
			for w_pos in sorted(chromosomes_positions_values[c]):
				if index%2==0:
					values.append(chromosomes_positions_values[c][w_pos][1])
				index+=1
	else:
		for c in chromosomes_positions_values:
			index=0
			for w_pos in sorted(chromosomes_positions_values[c]):
				if index%2==0:
					values.append(chromosomes_positions_values[c][w_pos])
				index+=1

	return values

###################################################################################
def check_depth_criteria_het(ad1, ad2):

	if ((ad1+ad2) >= 20) and (((ad1/(ad1+ad2)) >= 0.20) or ((ad1/(ad1+ad2)) <= 0.80)):
		return True
	else:
		return False

###################################################################################
def check_allelic_depth_criteria(ad1,ad2):

	''' Input: total depth for site'''
	#Filter SNPS based on total depth (>=20)'''
	if (ad1+ad2) >= 20 :
		return True
	else:
		return False

##################################################################################################
def filter_based_on_separator(genotype_information):

	#Check that separator is as expected and return list of genotypes
	if "/" in genotype_information[2][0]:
		no_alleles_for_site=len(genotype_information[2][0].split("/"))
		alleles_for_site=genotype_information[2][0].split("/")
		return alleles_for_site

	elif "|" in genotype_information[2][0]:
		no_alleles_for_site=len(genotype_information[2][0].split("|"))
		alleles_for_site=genotype_information[2][0].split("|")
		return alleles_for_site

	else:
		print(f"Position {genotype_information[1]} on chromosome {genotype_information[0]} does not have expected genotype format: GT={genotype_information[2][0]}!")
		exit(1)

############################################################################################################
def check_biallelic(alleles):

	if len(alleles) > 2 or len(alleles)<2:
		return False
	else:
		return True

################################################################################################################################################
def extract_snp_variants_for_sample_util(filename, no_snps_uncalled_vcf, genotype_info_chrom_positions, column_no, selected_chrom_to_use):

	'''
	reads vcf and stores SNPs.
	'''
	if selected_chrom_to_use==0:
		fh_vcf_in=open(filename, "r")
		for line in fh_vcf_in:
			line_striped = line.rstrip("\n")

			if not line_striped.startswith("#"):
				#total no snps in vcf for sample (including uncalled sites)
				no_snps_uncalled_vcf+=1
				elements_of_columns=line_striped.split("\t")
				# Assumes that genotype column of selected sample (i.e. column_no) starts with the following info GT:AD:DP (e.g. 0/1:12,6:18). Takes account only biallelic sites for plotting (for heterozygosity estimates)
				genotype_info=elements_of_columns[int(column_no)].split(":")
				#SNP no. as key. Assumes that chrom number is given on first column of VCF and chromosome position is given on second column of VCF. Value structure -> [chrom, position, [GT, AD, DP]]
				genotype_info_chrom_positions[no_snps_uncalled_vcf]=[elements_of_columns[0], elements_of_columns[1], genotype_info]
		fh_vcf_in.close()

	else:
		fh_vcf_in=open(filename, "r")
		for line in fh_vcf_in:
			line_striped = line.rstrip("\n")
			elements_of_columns=line_striped.split("\t")

			if not line_striped.startswith("#") and elements_of_columns[0] in selected_chrom_to_use:
				#total no snps in vcf for sample (including uncalled sites)
				no_snps_uncalled_vcf+=1
				# Assumes that genotype column of selected sample (i.e. column_no) starts with the following info GT:AD:DP (e.g. 0/1:12,6:18). Takes account only biallelic sites for plotting (for heterozygosity estimates)
				genotype_info=elements_of_columns[int(column_no)].split(":")
				#SNP no. as key. Assumes that chrom number is given on first column of VCF and chromosome position is given on second column of VCF. Value structure -> [chrom, position, [GT, AD, DP]]
				genotype_info_chrom_positions[no_snps_uncalled_vcf]=[elements_of_columns[0], elements_of_columns[1], genotype_info]
		fh_vcf_in.close()

	return (no_snps_uncalled_vcf, genotype_info_chrom_positions)
