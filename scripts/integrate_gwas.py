#!/usr/bin/python
###########################################################################################
# A script that parses the input from CSV file with GWAS data (published or unpublished)  #
# and SNPs that are in LD with those from GWAS                                            #
#                                               #
# This script takes 2 input parameters:                                                   #
#                                                               #
# (a) GWAS input file (should be TAB-delimited with columns with column names):           #
#                                                  #
# 1. PUBMED_id                                                                            #
# 2. Disease/Trait_name                                                                   #
# 3. Chr_id                                                                               #
# 4. Chr_position                                                                         # 
# 5. Mapped_gene                                                                          #
# 6. Strongest_SNP-risk_allele                                                            #
# 7. SNP name                                                                             #
# 8. P_value                                                                              # 
#                                                                                         #
# All columns need to exist, but if some information is missing/unavailable,              #
# values in these fields should be set to 0                                    #
#                                                              # 
# (b) Input file with LD information can be easily created following few steps using 
#    1) From your GWAS input file, Extract only a column with SNP names as a separate file (without column name)
#    2) Go to online tool HaploReg accessible through http://www.broadinstitute.org/mammals/haploreg/haploreg_v2.php to generate the list of LD SNPs 
#    3) Once the page has finished with data loading in your browser, click on your page (some blank part of it) and use CTRL+a to mark all the text
#    4) Copy the text using either CTRL+c or by Right-Click/Copy on the page
#    5) Save copied text in a separate file which will serve as an input
# and providing a TXT file using only SNP names from GWAS input file (without column names) #                 
#############################################################################################

import mySQLConnect
from mySQLConnect import mysql_connect
import sys
import csv
import re
import argparse
import os
import multiprocessing 
from os.path import basename

from load_gwas_data import load_gwas
from load_haploreg_data import load_haploreg
from remove_snps_with_pos0 import remove_incomplete_data
from filter_snps2 import filter_dhlSNPs
from count_lf import count_SNPs
from create_blocks_across_diseases import define_GRRs
from assign_snp_regions_to_bins2 import bin_GRRs

# Connect to GARLIC DB
con=mysql_connect()
cursor=con.cursor()

def integrateGWAS(*args):

    if args[0].gwas_input is None or args[0].ld_input is None:
        print "\nUsage: garlic integrateGWAS <path_to_input_gwas_file> <path_to_input_ld_file>"
        print "Error: To few arguments"
        sys.exit(1)
    else:
        gwas_path_input=args[0].gwas_input[0]
        ld_path_input=args[0].ld_input[0]

    print "\n### Reading GWAS SNP data ###\n"    
    load_gwas(gwas_path_input)

    print "\n### Reading LD SNP data ###\n"
    load_haploreg(ld_path_input)

    print "\n### Removing inconsistent (incomplete) data ###\n"
    remove_incomplete_data()

    print "\n### Adjusting SNP data across diseases ###\n"
    filter_dhlSNPs()

    print "\n### Counting the number of SNPs (leaders & followers) in each disease ###\n"  
    count_SNPs()

    print "\n### Defining GRRs ###\n"
    define_GRRs()

    print "\n### Assigning GRRs to bins ###\n"
    bin_GRRs()
