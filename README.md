
![Alt text](chromgenomeplot.svg)

## Description
ChromGenomePlot provides a set of functions for analysis and visualization of diverse genomic features along the chromosomes (or scaffolds) of reference genome assemblies. It is primarily designed for analyzing highly contiguous (few scaffolds) or chromosome-level reference genome assemblies. It leverages 
tools and functions from matplotlib, seaborn, pandas, numpy and other popular python libraries. It can plot and calculate overall percent heterozygosity, major allele frequency along chromosomes using windows of specified size (kb), major allele frequency spectra (MAF). It can also plot coverage depth, GC content and other feature content based on supplied CHROM coordinate files (i.e., feature coordinate files) along chromosomes/scaffolds using windows of specified length.

## Usage
```
git clone https://github.com/alexvasilikop/ChromGenomePlot.git
cd ChromGenomePlot/
python3 chromgenomeplot.py -h
```
The subdirectory ```analyses_scripts/``` and the ```chromgenomeplot.py``` script should be in the current working directory (where the analysis is run). The script looks for the presence of the analysis scripts subdirectory in the current working directory and throws an error if it does not find it.

Depending on the analysis the following input files are required: 1) assembly fasta, 2) sorted BAM file of reads mapped onto the assembly, 3) Results of variant calling in VCF format that includes only SNPs (it should be in the format produced by GATK), 4) Genome feature annotation in CHROM coordinate format (e.g., genes or repeats). 

The CHROM format requires 3 columns and is 1-based-closed set of feature coordinates i.e., ```[start-pos, end-position]```. This in contrast to BED file which is zero-based open-end ```[start, end)```. It is better that the coordinates are sorted according to their start position The chrom file can be extracted from a GFF3 file with identical types of features as follows by running:
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
The last column is optional. Each line represents the coordinates of different features of the same type (e.g., different genes)

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
