#!/usr/bin/env python3

from collections import defaultdict
import matplotlib.pyplot as plt
import argparse

class Coverage_File():

	def __init__(self, filename):

		self.filename = filename
		#Dictionary of tuples: chromosome -> (position, depth)
		self.chromosomes_and_depths = defaultdict(lambda: [])
		self.region_windows_depths = defaultdict(lambda: 0)

	def read_coverage(self, chromosome):

		depths_file=self.filename

		with open(depths_file, "r") as fh_in:

			for line in fh_in:

				if line.startswith(chromosome):
					line = line.rstrip("\n")
					elements = line.split(sep="\t")
					self.chromosomes_and_depths[elements[0]].append((int(elements[1]), int(elements[2])))

	def get_average_coverages_for_region(self, chromosome, start, end, bin_size):

		start_window = start
		depth_total_window = 0
		count_bases_window = 0
		target=0

		for p, d in self.chromosomes_and_depths[chromosome]:

			if p == start:
				target=1
				depth_total_window += d
				count_bases_window += 1

			elif target==1:
				depth_total_window += d
				count_bases_window += 1

				if p%bin_size==0:

					average_depth = int(depth_total_window/bin_size)
					self.region_windows_depths[start_window] = average_depth
					self.region_windows_depths[p] = average_depth
					depth_total_window = 0
					start_window= p+1
					count_bases_window = 0

					if p==end:
						target=0
						break

				elif p == end:
					average_depth = int(depth_total_window/count_bases_window)
					self.region_windows_depths[start_window] = average_depth
					self.region_windows_depths[p] = average_depth
					target=0
					break

	def  return_values_for_plotting(self):

		return self.region_windows_depths

####################################################################################################################""
class Plot():

	def __init__(self, fig_name):

		self.fig_name = fig_name

	def generate_plot(self, species, chromosome, chrom_depths):

		fig, axs = plt.subplots(1)
		fig.suptitle(chromosome+" - "+species, fontsize=10)
		tuples = sorted(chrom_depths.items()) # sorted by key, return a list of tuples
		x, y=zip(*tuples) # unpack a tuple into two list pairs
		plt.locator_params(axis='y', nbins=5)
		axs.grid()
		axs.set_ylabel('Depth of coverage')
		axs.set_xlabel('Position in selected chromosome region')
		axs.plot(x, y, color="steelblue")
		axs.fill_between(x=x, y1=y, color='blue', alpha=.1)
		fig.savefig(self.fig_name)
		plt.show()

########################################################################################################################################################################
def main():

	print("\n###############################################################################################################################################")
	print("Plotting coverage depth for a specified region on a specific chromosome of the reference genome...")
	print("##################################################################################################################################################\n")

	parser=argparse.ArgumentParser(description="Plot coverage along selected region on a specified chromosome. Requires depths.txt file as an argument.")
	parser.add_argument("in_depths_file", help="Depths file generated with e.g. SAMtools with three comumns tab delimited (chromosome position depth)")
	parser.add_argument("out_plot", help="Output plot")
	parser.add_argument("bin_size", help="Size of bins (bp) for coverage plot (average coverage plotted for each bin)")
	parser.add_argument("start", help="Start position of target region")
	parser.add_argument("end", help="End position of target region")
	parser.add_argument("chromosome", help="Target chromosome")
	parser.add_argument("species", help="Species name/library for title of the plot")
	parser.usage = 'python3 chromgenomeplot.py cov_depth_region [positional arguments]'

	args = parser.parse_args()

	#Assign bam/sam file as an instance of the Alignment Class
	my_depths=Coverage_File(filename=args.in_depths_file)

	#read BAM/SAM file and calculate per base coverages -> generates depths.txt file
	print("Reading depths file  ...\n")
	my_depths.read_coverage(args.chromosome)

	#calculates coverage per window along the chromosomes
	print("Getting average coverage depth for each window in specified region  ...\n")
	my_depths.get_average_coverages_for_region(args.chromosome,int(args.start),int(args.end),int(args.bin_size))

	#generate plot object
	my_plot=Plot(fig_name=args.out_plot)
	print("Plotting coverage depth for each window in specified region  ...\n")
	my_plot.generate_plot(args.species, args.chromosome, my_depths.return_values_for_plotting())

######################################################################################################################################################################
if __name__ == '__main__':
	main()
