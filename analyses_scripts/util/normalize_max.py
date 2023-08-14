#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def normalize(values_coverage_chromosome, max_plot, average_depth, no_avg_cov):

	if not no_avg_cov:

		# if --no_avg_cov is not activated then the coverage is normalized according to the average across all chromosomes and multiplied by max_snp or max_het (divided by 2).
		#For example if max_het option is 10, then the normalized coverage is multiplied by 10/2=5 for the plot)

		norm_average=[]
		for val in values_coverage_chromosome:
			norm_average.append((val/average_depth)*(max_plot/2))

		return norm_average

	else:

		#if --no_avg_cov is activated then the coverage is not normalized according to the average across all chromosomes and only based on the max_snp or max_het value.
		#This could be a problem if on some chromosomes you have regions with outlier coverage depth (e.g., due to repetitive elements)

		norm_max=[]
		for val in values_coverage_chromosome:
			norm_max.append(((val-min(values_coverage_chromosome))/(max(values_coverage_chromosome)-min(values_coverage_chromosome)))*max_plot)

		return norm_max

