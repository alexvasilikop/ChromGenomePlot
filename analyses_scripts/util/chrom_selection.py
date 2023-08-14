#!/usr/bin/env python3

def chromosome_selection(argument):

	#Check if all or some chromosomes will be used for the output
	selected_chrom_to_use=0
	if not (argument):
		print("All chromosomes will be used for calculations and the output plot...\n")
	else:
		selected_chrom_to_use=argument.split(",")

		print("Selected chromosomes:")
		for i in selected_chrom_to_use:
			print(i)
		print("\n"+200*"-")

	return selected_chrom_to_use