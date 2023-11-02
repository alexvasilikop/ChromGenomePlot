#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from util import chrom_selection, feature_utils
from feature import Assembly_FASTA, Plot
from collections import defaultdict
import matplotlib.pyplot as plt
import argparse

###########################################################################################################################################################################################################################
# Plots GC content and feature content along chromosomes using bed files of annotations (each bed annotation must contain sorted annotations by genome position annotations)
##############################################################################################################################################################################################################

class Assembly_FASTA(Assembly_FASTA):

	def calculate_average_bed_values_chromosomes(self, bin_size):
		
		for chromosome in self.scaffolds_seqs.keys():

			print("Working on chromosome: "+chromosome+" ...")
			#Coverage-specific values
			start_window = 1
			count_bases_window = 0
			count_AT=0
			count_GC=0
			length_chromosome = len(self.scaffolds_seqs[chromosome])
			position=0

			for nuc in self.scaffolds_seqs[chromosome]:
				count_bases_window += 1
				if nuc.upper() in "AT":
					count_AT+=1
				elif nuc.upper() in "GC":
					count_GC+=1

				if (position+1)%bin_size==0:

					#Add depth, snp percent of start position
					self.chromosomes_windows_coverage_numbers[chromosome][start_window] = ((count_GC/(count_GC+count_AT))*100)
					#Add depth , snp percent  of end position
					self.chromosomes_windows_coverage_numbers[chromosome][position+1] = ((count_GC/(count_GC+count_AT))*100)
					start_window= (position+1)+1
					count_bases_window = 0
					count_AT=0
					count_GC=0

				#if reaching the end of chromosome before the windowstep is completed
				elif position == length_chromosome-1:

					#Add depth, snp percent of start position
					self.chromosomes_windows_coverage_numbers[chromosome][start_window] = ((count_GC/(count_GC+count_AT))*100)
					#Add depth , snp percent  of end position
					self.chromosomes_windows_coverage_numbers[chromosome][position+1] = ((count_GC/(count_GC+count_AT))*100)

				position+=1

###############################################################################################################################################################################
class Plot(Plot):

	def plot_features_chromosomes(self, features_chroms, max_len, species, no_fill):

		plt.rcParams["font.family"]= "Arial"

		fig, axs = plt.subplots(nrows=len(features_chroms), ncols=1, figsize=(14,10), sharey=True, sharex=True, tight_layout=True)
		title = fig.suptitle('GC content along chromosomes'+" - "+species, fontsize=10)
		x_title = fig.supxlabel("Position on chromosome (Mb)")
		y_title = fig.supylabel("Percentage of bases covered in chromosome window (%)")
		no=0

		if len(features_chroms.keys())>1:
			for i in features_chroms.keys():

				#Unpack dictionary -> creates lists of values for plotting
				position = []
				y_feat2 = []

				for p, value in sorted(features_chroms[i].items()):
					position.append(p)
					y_feat2.append(value)

				l2 = axs[no].plot(position, y_feat2, color="orange", label="GC content (%)")
				if not no_fill:
					axs[no].fill_between(x=position, y1=y_feat2, color='orange', alpha=.1)

				axs[no].locator_params(axis='y', nbins=4)
				axs[no].set_title(i, fontsize=8)
				axs[no].grid()
				no+=1

		else:
			position = []
			y_feat1_coverage = []
			y_feat2 = []

			for p, value in features_chroms[list(features_chroms.keys())[0]].items():
				position.append(p)
				y_feat2.append(value)
			
			l2 = axs.plot(position, y_feat2, color="orange", label="GC content (%)")
			if not no_fill:
				axs.fill_between(x=position, y1=y_feat2, color='orange', alpha=.1)
			axs.locator_params(axis='y', nbins=4)
			
			axs.set_title(list(features_chroms.keys())[0], fontsize=8)
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

		if len(features_chroms.keys())>1:
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
	print("Estimating GC content for chromosome windows and plotting along the chromosomes of the reference genome...")
	print("############################################################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot GC content along the chromosomes")
	parser.add_argument("in_assembly", help="Genome assembly fasta used to make the annotations.")
	parser.add_argument("out_plot", help="Output plot with GC content and feature content along the chromosomes")
	parser.add_argument("bin_size", help="Size of bins (bp) for plot (GC percentage plotted for each bin)")
	parser.add_argument("species", help="Species name/sequencing library for title of the plot")
	parser.add_argument('-sel_chrom', action='store', help="List of selected chromosomes to use (separated by \",\" without spaces, e.g.: \"chrom_1,chrom_2\" )")
	parser.add_argument('--no_fill', action='store_true', help="Use this flag to indicate no color filling between the lineplot and the x axis")

	parser.usage = 'python3 chromgenomeplot.py GC_content [positional arguments]'

	#Parsing args
	args = parser.parse_args()
	selected_chromosomes=chrom_selection.chromosome_selection(args.sel_chrom)

	#Parsing assembly fasta
	print("\nProcessing Assembly fasta file...")
	my_assembly = Assembly_FASTA(filename=args.in_assembly)
	print("Reading Assembly fasta file...")
	my_assembly.read_fasta(selected_chromosomes)
	print("Finished reading Assembly fasta file...\n"+200*"-"+"\n")

	#Calculates GC content (percent bases covered) per window along the chromosomes
	print("Extracting average GC content for windows of "+str(args.bin_size)+" bp ...")
	my_assembly.calculate_average_bed_values_chromosomes(bin_size=int(args.bin_size))
	print("Finished processing all files!\n"+200*"-"+"\n")
	
	#Return values for plotting
	print("Preparing the plot...")
	my_plot = Plot(fig_name=args.out_plot)
	my_plot.plot_features_chromosomes(features_chroms=my_assembly.chromosomes_windows_coverage_numbers, \
									  max_len=max([len(seq) for seq in my_assembly.scaffolds_seqs.values()]), \
									  species=args.species, \
									  no_fill=args.no_fill)
	print("All done!")

########################################################################################################################################################################
if __name__ == '__main__':
	main()

	


