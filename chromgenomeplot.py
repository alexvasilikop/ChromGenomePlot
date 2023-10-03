#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv
import subprocess
import os

##################################################################################################
def starting_message():

    print(2*"\n"+"######################################################################################")
    print("#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#")
    print("#<                                                                                  >#")
    print("#<     ChromGenomePlot: analysis and visualization of diverse genomic features      >#")
    print("#<             along the chromosomes of reference genome assemblies                 >#")
    print("#<                                                                                  >#")
    print("#<                        email: alexvasilikop@gmail.com                            >#")
    print("#<                                                                                  >#")
    print("#⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄⌄#")
    print("######################################################################################"+2*"\n")

##########################################################################################################"""
def check_scripts_dir():

    #check analyses dir exists in directory
    if os.path.isdir("analyses_scripts"):
        pass
    else:
        print("Directory \"analyses_scripts\" is not in the current working directory!\n")

######################################################################################################################################
def check_analysis_script(analysis_type):

    #check analyses script exists in "analyses_scripts" directory
    if os.path.isfile("analyses_scripts/"+analysis_type+".py"):
        pass
    else:
        print("\nScript "+analysis_type+".py"+" not in \"analyses_scripts\" directory!\n")
        exit(1)

################################################################################################################################
def help_individual_analyses(analysis_type):

    starting_message()
    #print type of analysis
    print("Analysis type:\n-> "+analysis_type+"\n")

    #check script is in the analysis_scripts directory
    check_analysis_script(analysis_type)

    #run help of invividual analysis
    command = "python3 "+"analyses_scripts/"+argv[1]+".py -h"
    proc = subprocess.Popen(command, shell=True)
    proc.communicate()

##################################################################################################################################
def main():

    analyses_types = ("cov_hist", "cov_depth", "cov_depth_region", "feature", "feat_2_feat", \
                      "cov_het", "cov_snps", "cov_snps_number", "snps_number", "snps_percent", \
                      "snps_het", "allele_freq_chrom", "allele_freq_dist", "allele_freq_d_sel", \
                      "GC_2_feat", "GC_content", "cov_vs_GC")

    #check analyses scripts directory
    check_scripts_dir()

    if len(argv)<2:
        starting_message()
        print("Incorrect number of arguments provided!\nGeneral usage:\npython3 chromgenomeplot.py [analysis_type] [positional_arguments]\n")
        print("Help message:\npython3 chromgenomeplot.py -h")

    elif "-h" in argv and len(argv)==2:

            #print help message (all different types of analyses)
            starting_message()
            print("General usage:\npython3 chromgenomeplot.py [analysis_type] [positional_arguments]\n")
            print("\nTo see options for each type of analysis:\npython3 chromgenomeplot.py [analysis_type] -h"+"\n\n")
            print("\t Types of analyses supported:\n\n")
            print(f"\t {analyses_types[0]:<20} Plot histogram of coverage depth for the reference genome assembly\n")
            print(f"\t {analyses_types[1]:<20} Plot coverage depth along the chromosomes of the reference genome\n")
            print(f"\t {analyses_types[2]:<20} Plot coverage depth of a specific region on a specified chromosome\n")
            print(f"\t {analyses_types[3]:<20} Plots percent coverage of one genomic feature in each window along chromosomes of the reference genome\n")
            print(f"\t {analyses_types[4]:<20} Plots percent coverage of two genomic features in each window along chromosomes of the reference genome\n")
            print(f"\t {analyses_types[5]:<20} Infer average coverage depth and SNP-based heterozygosity (percentage of heterozygous SNPs) and plot along chromosomes\n")
            print(f"\t {analyses_types[6]:<20} Plot SNP content (percent bases, %) and average depth along the chromosomes of the reference genome\n")
            print(f"\t {analyses_types[7]:<20} Plot SNP number and average coverage depth along the chromosomes of the reference genome\n")
            print(f"\t {analyses_types[8]:<20} Plot SNP number by chromosome window along the chromosomes of the reference genome\n")
            print(f"\t {analyses_types[9]:<20} Plot SNP percentage (% bases covered) by chromosome window along the chromosomes of the reference genome\n")
            print(f"\t {analyses_types[10]:<20} Infer SNP-based heterozygosity (% of heterozygous SNPs) and plot along the chromosomes of the reference genome\n")
            print(f"\t {analyses_types[11]:<20} Infer dominant allele frequency for chromosome windows and plot it along the chromosomes of the reference genome\n")
            print(f"\t {analyses_types[12]:<20} Make histogram (distribution plot) of dominant allele frequencies using a set of SNPs from a VCF file\n")
            print(f"\t {analyses_types[13]:<20} Make histogram (distribution plot) of dominant allele frequencies for a specified region of a selected chromosome\n")
            print(f"\t {analyses_types[14]:<20} Plot GC vs feature content along the chromosomes of the reference genome\n")
            print(f"\t {analyses_types[15]:<20} Plot GC content content along the chromosomes of the reference genome\n")
            print(f"\t {analyses_types[16]:<20} Plot GC content versus coverage depth of the chromosomes in the reference genome\n")
            
    else:

        if argv[1] not in analyses_types:

            #wrong type of analysis specified
            starting_message()
            print("Analysis \""+argv[1] + "\" not found!\n")
            print("Help message to show types of analyses:\npython3 chromgenomeplot.py -h\n")
            print("General usage:\npython3 chromgenomeplot.py [analysis_type] [positional_arguments]\n")
            print("To see options for each analysis:\npython3 chromgenomeplot.py [analysis_type] -h"+"\n")
            exit(1)

        else:

            #print help for individual analyses
            if "-h" in argv:
                help_individual_analyses(argv[1])

            else:
                starting_message()
                print("Analysis type:\n-> "+argv[1])

                #parse arguments provided to the analysis type
                arguments = [arg for arg in argv[2:]]
                
                if len(arguments)==0:
                    print("\n"+"No command line arguments were provided!")
                    print("\n"+"Usage for analysis type "+argv[1]+":\n"+"python3 chromgenomeplot.py"+" "+argv[1]+" [positional_arguments]\n")

                    #print general usage
                    print("General usage of chromgenomeplot:\npython3 chromgenomeplot.py [analysis_type] [positional_arguments]\n")

                    #help message
                    print("Help message to show all types of analyses:\npython3 chromgenomeplot.py -h\n")
                    print("To see options for "+argv[1]+" analysis:\npython3 chromgenomeplot.py "+argv[1]+" -h"+"\n")
                    exit(1)

                print("\n"+"The following command line arguments were provided:")
                for arg in argv[2:]:
                    print(arg)
                print("\n"+str(200*"-"))

                command = "python3 "+"analyses_scripts/"+argv[1]+".py "+" ".join(argv[2:])
                print("Running CMD: "+command+"\n")
                proc = subprocess.Popen(command, shell=True)
                proc.communicate()

#####################################################################################################################################################
if __name__ == '__main__':
    main()




   
