# ChromGenomePlot
ChromGenomePlot provides a set of functions for analysis and visualization of diverse genomic features along the chromosomes of reference genome assemblies. It is primarily designed for analyzing highly contiguous (few scaffolds) or chromosome-level genome assemblies. It leverages 
tools and functions from matplotlib, seaborn, numpy and other popular python libraries. It can plot and calculate overall percent heterozygosity and/or along chromosomes using windows of specified size (kb). It can also plot coverage depth, GC content and other features based on .bed files along the chromosome windows.

## Usage
```
python3 chromgenomeplot.py -h
```
The subdirectory ```analyses_scripts/``` should be copied in the same dir as chromgenomeplot.py
Depending on the analysis the following input files are required: 1) assembly fasta, 2) sorted BAM file of reads mapped onto the assembly, 3) Results of variant calling in VCF format that includes only SNPs, 4) Genome annotation in BED format. 


## Dependencies
coming soon..

## Description of different functionalities
coming soon..
