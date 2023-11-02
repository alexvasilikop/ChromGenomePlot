#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import shutil

################################################################################################
def run_samtools(file, depths_file) :

	#check SAMtools is on path
	if shutil.which("samtools"):
		print("SAMtools is installed and on path...")
	else:
		print("SAMtools is not installed or is not on path...")
		exit(1)
	
	cmd="samtools depth -aa"+" "+file
	fh_out = open(depths_file, "w")
	proc = subprocess.Popen(cmd, shell=True, stdout=fh_out)
	proc.communicate()
	fh_out.close()

##############################################################################################
def read_depths_file_util(depths_file, depth_total, selected_chrom_to_use, chromosomes_and_depths, genome_size):
	#all chromosomes
	#mthis method is faster than using pandas to read in a dataframe
	if selected_chrom_to_use==0:
		with open(depths_file, "r") as fh_in:
			for line in fh_in:
				line = line.rstrip("\n")
				elements = line.split(sep="\t")
				chromosomes_and_depths[elements[0]].append((int(elements[1]), int(elements[2])))
				depth_total+=int(elements[2])
				genome_size+=1

	#selection of chromosomes
	else:
		with open(depths_file, "r") as fh_in:

			for line in fh_in:
				line = line.rstrip("\n")
				elements=line.split("\t")

				if elements[0] in selected_chrom_to_use:
					chromosomes_and_depths[elements[0]].append((int(elements[1]), int(elements[2])))
					depth_total+=int(elements[2])
					genome_size+=1

	return (depth_total, chromosomes_and_depths, genome_size)
