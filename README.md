
![Alt text](https://github.com/alexvasilikop/ChromGenomePlot/blob/main/chromgenomeplot.svg)

## Description
ChromGenomePlot provides a set of functions for analysis and visualization of diverse genomic features along the chromosomes of reference genome assemblies. It is primarily designed for analyzing highly contiguous (few scaffolds) or chromosome-level reference genome assemblies. It leverages 
tools and functions from matplotlib, seaborn, pandas and other popular python libraries. It can plot and calculate overall percent heterozygosity, major allele frequency along chromosomes using windows of specified size (kb), major/minor allele frequency spectra (MAF). It can also plot coverage depth, GC content and other feature content based on .bed files along the chromosomes based on chromosome windows of specified length.

## Usage
```
python3 chromgenomeplot.py -h
```
The subdirectory ```analyses_scripts/``` should be copied in the same dir as ```chromgenomeplot.py``` or in the working directory.

Depending on the analysis the following input files are required: 1) assembly fasta, 2) sorted BAM file of reads mapped onto the assembly, 3) Results of variant calling in VCF format that includes only SNPs (it should be in the format produced by GATK), 4) Genome feature annotation in BED format (e.g., genes or repeats). 


## Dependencies
The tool requires installation of python 3. Installation of the following python libraries are required
```
matplotlib
seaborn
pandas
```

In addition ```Samtools``` should be installed and on path. To control the versions of libraries and software, it is recommended to install them within a conda envronment. It has been tested using the following versions:
```
matplotlib                3.5.2
seaborn                   0.11.2
numpy                     1.23.1
pandas                    1.4.3
samtools                  1.10 (using htslib 1.10.2)
```

## Description of different functionalities
coming soon..
