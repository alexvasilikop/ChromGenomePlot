#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
from Bio import SeqIO
import subprocess
from collections import defaultdict
import shutil

def make_BED_intervals(scaffolds_seqs, bin_size, outfile):
	#make intervals in BED format for the chromosome windows

	fh_out_intervals=open(outfile, "w")
	
	for chrom in scaffolds_seqs:
	    l_chrom=len(scaffolds_seqs[chrom])

	    for i in range(l_chrom):

	    	if i%bin_size==0:

	    		if (i+1+bin_size)>l_chrom:
	    			fh_out_intervals.write(chrom+"\t"+str(i)+"\t"+str(l_chrom)+"\n")
	    		else:
	    			fh_out_intervals.write(chrom+"\t"+str(i)+"\t"+str(i+bin_size)+"\n")
	fh_out_intervals.close()

######################################################################################################
def CHROM_read_util(filename, chrom_features, selected_chrom_to_use):
	
	#all chromosomes read chrom coordinates
	if selected_chrom_to_use==0:
		with open(filename, "r") as fh_in:

			for line in fh_in:
				line = line.rstrip("\n")
				elements = line.split("\t")
				#assumes 1-based, end-closed [start, end] CHROM-formatted data.
				chrom_features[elements[0]].extend([int(elements[1]), int(elements[2])])

	#selection
	else:
		with open(filename, "r") as fh_in:

			for line in fh_in:
				line = line.rstrip("\n")
				elements = line.split("\t")

				if elements[0] in selected_chrom_to_use:
					#assumes 1-based, end-closed [start, end] CHROM-formatted data.
					chrom_features[elements[0]].extend([int(elements[1]), int(elements[2])])

	return chrom_features

#########################################################################################
def fasta_read_util(filename, scaffolds_seqs, selected_chrom_to_use):
	#read fasta
	if selected_chrom_to_use==0:
		with open(filename) as fh_in:
			for record in SeqIO.parse(fh_in, "fasta"):
				scaffolds_seqs[record.id]=record.seq

	#keep only selected chromosomes
	else:
		with open(filename) as fh_in:
			for record in SeqIO.parse(fh_in, "fasta"):
				scaffolds_seqs[record.id]=record.seq

		#keep only selected chromosomes
		try:
			IDs=list(scaffolds_seqs.keys())
			for i in IDs:
				if i not in selected_chrom_to_use:
					del scaffolds_seqs[i]
		except (KeyError):
			print("Key {i} Not In FASTA")

	return scaffolds_seqs

########################################################################################################
def total_window_numbers_or_coverage(bin_size, scaffolds_seqs, feature, label, mode_numbers):
	#relies on bedtools to estimate numbers of feature in chromosome intervals
	#Writes intermediate files (interval bed and feature bed)
	print("\n\n"+100*"#"+"\n# Starting coverage calculations with BEDTools for "+label+" ###\n"+100*"#")

	#check BEDTools is on path
	if shutil.which("bedtools"):
		print("BEDTools is installed and on path...")
	else:
		print("BEDTools is not installed or is not on path...")
		exit(1)

	#Interval BED write
	print("\n## Constructing intervals in BED format for BEDTools ("+label+") ... ##\n")
	make_BED_intervals(scaffolds_seqs, bin_size, outfile="bed_intervals_"+label+".txt")

	#Feature BED write
	print("## Constructing BED feature interval file for BEDTools ("+label+") ... ## \n")
	fh_out_feature_intervals=open("feature_intervals_"+label+".txt", "w")
	for chrom in feature:
		index=0

		while index < len(feature[chrom])-1:
			fh_out_feature_intervals.write(chrom+"\t"+str((feature[chrom][index]))+"\t"+str((feature[chrom][index+1]))+"\n")
			index+=2
	fh_out_feature_intervals.close()

	#Run bedtools and write output in file
	print("## Running BEDTools coverage and write output in file ("+label+") ... ##\n")
	cmd=("bedtools coverage -a bed_intervals_"+label+".txt -b feature_intervals_"+label+".txt")
	bedtools_out=open("bedtools_coverage_and_numbers_"+label+".txt", "w")
	process = subprocess.Popen(cmd, shell=True, stdout=bedtools_out)
	process.communicate()
	bedtools_out.close()

	#Read file and return numbers or coverage of features for intervals
	scaffolds_and_feat_numbers= defaultdict(lambda: {})

	with open("bedtools_coverage_and_numbers_"+label+".txt", "r") as fh_in:

		for line in fh_in:
			line = line.strip()
			elements=line.split("\t")

			#create data structure for plotting
			if mode_numbers:
				scaffolds_and_feat_numbers[elements[0]][int(elements[1])+1]=int(elements[3])
				scaffolds_and_feat_numbers[elements[0]][int(elements[2])]=int(elements[3])

			else:
				scaffolds_and_feat_numbers[elements[0]][int(elements[1])+1]=float(elements[6])*100
				scaffolds_and_feat_numbers[elements[0]][int(elements[2])]=float(elements[6])*100

	print("# Finished calculations with BEDTools for "+label+" ###\n"+100*"#")
	return scaffolds_and_feat_numbers

########################################################################################################
def coverage_or_numbers_per_window(chroms_and_feat_content, mode_numbers):
	#returns average content (coverage % or number of features per window)
	list_values=[]

	if mode_numbers:
		for chrom in chroms_and_feat_content.keys():
			index=0
			for position in chroms_and_feat_content[chrom]:
				if index%2==0:
					list_values.append(int(chroms_and_feat_content[chrom][position]))

	else:
		for chrom in chroms_and_feat_content.keys():
			index=0
			for position in chroms_and_feat_content[chrom]:
				if index%2==0:
					list_values.append(int(chroms_and_feat_content[chrom][position]))
	return list_values

#########################################################################################################
def get_GC_content_for_windows(chroms_and_feat_content, scaffolds_seqs, bin_size):

	#make new data structure for adding GC percent as a second element in the value of the dictionary (similar to second feature)
	new_features=defaultdict(lambda: {})
	
	for c in chroms_and_feat_content:
		index=0
		for pos in sorted(chroms_and_feat_content[c]):
			
			if index==0:
				G_content=scaffolds_seqs[c][int(pos-1):int(pos-1+bin_size)].upper().count('G')
				C_content=scaffolds_seqs[c][int(pos-1):int(pos-1+bin_size)].upper().count('C')
				GC_percent=((G_content+C_content)/bin_size)*100

				new_features[c][pos]=[chroms_and_feat_content[c][pos], GC_percent]
				new_features[c][pos-1+bin_size]=[chroms_and_feat_content[c][pos-1+bin_size], GC_percent]

			elif index>=1 and index%2==0 and index<(len(chroms_and_feat_content[c])-2):
				G_content=scaffolds_seqs[c][int(pos-1):int(pos-1+bin_size)].upper().count('G')
				C_content=scaffolds_seqs[c][int(pos-1):int(pos-1+bin_size)].upper().count('C')
				GC_percent=((G_content+C_content)/bin_size)*100

				new_features[c][pos]=[chroms_and_feat_content[c][pos], GC_percent]
				new_features[c][pos-1+bin_size]=[chroms_and_feat_content[c][pos-1+bin_size], GC_percent]

			elif index==(len(chroms_and_feat_content[c])-2):
				fragment_length=len(scaffolds_seqs[c][int(pos-1):len(scaffolds_seqs[c])])
				G_content=scaffolds_seqs[c][int(pos-1):int(pos-1+fragment_length)].upper().count('G')
				C_content=scaffolds_seqs[c][int(pos-1):int(pos-1+fragment_length)].upper().count('C')
				GC_percent=((G_content+C_content)/fragment_length)*100

				new_features[c][pos]=[chroms_and_feat_content[c][pos], GC_percent]
				new_features[c][pos-1+fragment_length]=[chroms_and_feat_content[c][pos-1+fragment_length], GC_percent]
			index+=1
	
	return new_features
