#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Generates coverage depth histogram plot (per base coverage depth histogram) from an input BAM file (sorted).
#It requires that samtools is in the path. If depths.txt file is already in the directory this is used by default.

import subprocess
import argparse
import matplotlib.pyplot as plt
from util.samtools_depths_util import run_samtools
import matplotlib.ticker as mticker
import pandas as pd
import seaborn as sns
import os

class Alignment():

	def __init__(self, filename):

		self.filename = filename
		self.depths_df = pd.DataFrame()

	def read_depths_file(self, depths_file):

		self.depths_df = pd.read_csv(depths_file, sep="\t", header = None)
		self.depths_df.columns = ['chromosome', 'position', 'depth']

	def read_alignment(self):

		depths_file="depths.txt"

		if os.path.isfile(depths_file):
			print("Found existing depths file! Samtools depth function will be skipped..")
			print("Using previously inferred depths.txt file..")
			self.read_depths_file(depths_file)

		else:
			depths_file="depths.txt"
			#samtools must be on current working directory
			print("Running samtools depth...")
			run_samtools(self.filename, depths_file)
			self.read_depths_file(depths_file)

	def plot_histogram_coverage(self, out_plot, no_bins, max_depth):

		sns.set_style("darkgrid")

		#Set upper value limit for histogram
		sns.histplot(data=self.depths_df[self.depths_df['depth'] <= int(max_depth)]['depth'], bins=int(no_bins), kde = False, color='steelblue')

		#set titles
		x_title = plt.xlabel('Coverage depth')
		y_title = plt.ylabel('No. of bases')
		title = plt.title('Coverage plot')

		current_values = plt.gca().get_yticks().tolist()
		plt.gca().yaxis.set_major_locator(mticker.FixedLocator(current_values))
		plt.gca().set_yticklabels(['{:.0f}'.format(x) for x in current_values])
		plt.grid(True)
		plt.tight_layout()
		plt.savefig(out_plot, bbox_extra_artists=(title, x_title, y_title))
		plt.show()

##############################################################################################################
def main():

	print("\n############################################################################################################################################################")
	print("Inferring histogram of base coverage depths for the reference genome...")
	print("############################################################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Calculate coverage depth histogram from BAM file. Requires samtools on path")
	parser.add_argument("in_alignment", help="Genome alignment file in BAM format")
	parser.add_argument("out_plot", help="Output file in png format to save the plot")
	parser.add_argument("no_bins", help="Number of bins to use for the histogram plot")
	parser.add_argument("max_depth", help="Maximum depth to plot in the histogram")

	parser.usage="python3 chromgenomeplot.py cov_hist [positional_arguments]"
	args = parser.parse_args()

	#Read bam/sam file
	my_alignment=Alignment(filename=args.in_alignment)

	#read BAM/SAM file and calculate per base coverages
	print("Reading input BAM file and generating depths.txt file...")
	my_alignment.read_alignment()
	
	#plot histogram of depths
	print("Generating plot...")
	my_alignment.plot_histogram_coverage(args.out_plot, args.no_bins, args.max_depth)
	print("All done!")

######################################################################################################################
if __name__ == '__main__':
	main()

	
