#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
from Bio import SeqIO

######################################################################################################
def bed_read_util(filename, chrom_features, selected_chrom_to_use):

	#all chromosomes
	if selected_chrom_to_use==0:
		with open(filename, "r") as fh_in:

			for line in fh_in:
				line = line.rstrip("\n")
				elements = line.split("\t")
				#assumes 0-based, half-open [start-1, end) extended BED-formatted data.
				positions = [x for x in range(int(elements[1]), int(elements[2]))]
				chrom_features[elements[0]].extend(positions)

	#selection
	else:
		with open(filename, "r") as fh_in:

			for line in fh_in:
				line = line.rstrip("\n")
				elements = line.split("\t")

				if elements[0] in selected_chrom_to_use:
					#assumes 0-based, half-open [start-1, end) extended BED-formatted data.
					positions = [x for x in range(int(elements[1]), int(elements[2]))]
					chrom_features[elements[0]].extend(positions)

	return chrom_features

###################################################################################################
def fasta_read_util(filename, scaffolds_seqs, selected_chrom_to_use):

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


#########################################################################################
def get_max_length(scaffolds_seqs):

	max_len = 0
	for seq in scaffolds_seqs.values():

		if len(seq)>max_len:
			max_len=len(seq)

	return max_len

#############################################################################################""
def total_window(start_window, p, chromosome, feature_pos):

	#count feature for window
	no_window=0

	for i in feature_pos:
		if (i>=start_window) and (i<=p):
			no_window+=1

	return int(no_window)
