
![Alt text](chromgenomeplot.svg)

## Description
ChromGenomePlot is a command-line tool written in python 3 that provides a set of functions for analysis and visualization of diverse genomic features along the chromosomes (or scaffolds) of reference genome assemblies. It is primarily designed to analye highly contiguous (few scaffolds) or chromosome-level reference genome assemblies. It leverages 
tools and functions from matplotlib, seaborn, pandas, numpy and other popular python libraries. It can plot and calculate overall percent heterozygosity, major allele frequency along chromosomes using windows of specified size (kb) as well as major allele frequency spectra (MAF). It can also plot coverage depth, GC content and other feature content based on supplied CHROM coordinate files (i.e., feature coordinate files) along chromosomes/scaffolds using chromosome windows of specified length.

## Usage
```
git clone https://github.com/alexvasilikop/ChromGenomePlot.git
cd ChromGenomePlot/
python3 chromgenomeplot.py -h
```
This will return all currently supported types of analysis. To print help message for each specific type of analyses type (some examples):
```
python3 chromgenomeplot.py cov_hist -h
python3 chromgenomeplot.py cov_het -h
python3 chromgenomeplot.py feature -h
```
This will return the help message and required positional (ordered arguments) for each analysis type. For example:
```
python3 chromgenomeplot.py feature -h
```
will return:
```

######################################################################################
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#<                                                                                  >#
#<     ChromGenomePlot: analysis and visualization of diverse genomic features      >#
#<             along the chromosomes of reference genome assemblies                 >#
#<                                                                                  >#
#<                        email: alexvasilikop@gmail.com                            >#
#<                                                                                  >#
#⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄#
######################################################################################


Analysis type:
-> feature


############################################################################################################################################################
Estimating feature content (one feature) for chromosome windows and plotting content along the chromosomes of the reference genome...
############################################################################################################################################################

usage: python3 chromgenomeplot.py feature [positional arguments]

Plot content of one feature (e.g., genes, coding sequences or repeats) along the chromosomes using a CHROM file of genome annotations

positional arguments:
  in_assembly           Genome assembly fasta used to make the annotations
  feature_1             CHROM coordinate file with annotated feature on the reference genome
  feature_1_name        Name of feature for plot label without spaces (e.g., coding_sequences, introns)
  out_plot              Output plot with feature content along chromosomes.
  bin_size              Size of bins (bp) for plot (percent bases covered by feature plotted for each bin)/or number of features if --numbers
  species               Species name / sequencing library for title of the plot

optional arguments:
  -h, --help            show this help message and exit
  -sel_chrom SEL_CHROM  List of selected chromosomes to use (separated by "," without spaces, e.g.: "chrom_1,chrom_2" )
  --no_fill             Use this flag to indicate no color filling between the lineplot and the x axis
  --numbers             Use this flag to indicate numbers of features instead of coverage of bases in chromosome windows
```

The subdirectory ```analyses_scripts/``` and the ```chromgenomeplot.py``` script should be in the current working directory (where the analysis is run). The script looks for the presence of the ```analysis_scripts``` subdirectory in the current working directory and throws an error if it does not find it. The tool has only been tested on Linux (Ubuntu v. 20)

## Input files
Depending on the analysis type, the following input files may be required: 1) assembly fasta, 2) sorted BAM file of reads mapped onto the assembly, 3) Results of variant calling in VCF format that includes only SNPs (it should be in the format produced by GATK4, see below for more details), 4) Genome feature annotation in CHROM coordinate format (e.g., genes or repeats). 

The CHROM format (required for feature-related functions) consists of three columns and is 1-based closed-end set of feature coordinates i.e., ```[start-pos, end-position]```. This in contrast to BED file which is zero-based open-end ```[start, end)```. It is better that the coordinates are sorted according to their start position. The chrom file can be extracted from a GFF3 file as follows by running:
```
cut -f1,3,4,5 chrom.gff3 |while read scaf type start end;do echo -e "$scaf\t$start\t$end\t$type";done |grep "gene"

Scaffold_11	16	444	gene
Scaffold_11	696	1259	gene
Scaffold_11	973	2020	gene
Scaffold_11	1996	2778	gene
Scaffold_11	2971	4033	gene
Scaffold_11	4945	7143	gene
Scaffold_11	7270	8945	gene
Scaffold_11	9004	9861	gene
Scaffold_11	9864	11177	gene
```
The last column is optional. Each row represents the coordinates of different features of the same type (e.g., different genes)

## Dependencies
The tool requires installation of python 3. Installation of the following python libraries are required
```
matplotlib
seaborn
pandas
numpy
Biopython
```

In addition ```SAMtools``` and/or ```BEDTools```should be installed and on path (depending on which analyses are performed). To control the versions of libraries and software, it is recommended to install them within a conda envronment. The pipeline has been tested using the following software and package versions:
```
matplotlib                3.5.2
seaborn                   0.11.2
numpy                     1.23.1
pandas                    1.4.3
SAMtools                  1.10 (using htslib 1.10.2)
BEDTools                  2.27.1
Biopython                 1.78
```

## Description of different functionalities
coming soon..
