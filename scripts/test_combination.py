#!/usr/bin/python
###################################################
#   A MAIN SCRIPT THAT STARTS ALL OTHER SCRIPTS   #
###################################################

import mySQLConnect
from mySQLConnect import mysql_connect
import sys
import csv
import re
import argparse
import os
import multiprocessing 
from os.path import basename
from load_new_map import load_map

from pcount_hits2 import main_pcount_hits    
from pcalculate_score4 import main_pcalculate_score    
from calculate_significance_light import main_calculate_significance    

con=mysql_connect()
cursor=con.cursor()

# Parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", nargs='+', type=str, help="Path to the input file (if not in DB)")
parser.add_argument("-s", nargs='+', type=int, help="Number of iterations for simulated score")
parser.add_argument("-d", nargs='+', type=str, help="Input disease id")
#parser.add_argument("-m", nargs='+', type=str, help="Regulatory map id (if already in DB)")

args=parser.parse_args()

if args.i is None or args.d is None:
#print "\nUsage: python start_test_combinations.py -i [<path_to_the_enhancer_map>|<map_id>] -d <disease id> [-s <num_iterations>]"
    print "\nUsage: python start_test_combinations.py -i <path_to_the_enhancer_map> -d <disease id> [-s <num_iterations>]"
    print "\nOptional parameters:"
    print "\t-s\tNumber of iterations for simulated score. Default value is 100000."
    print "Error: To few arguments"
    sys.exit(1)
else:
    if (args.i is not None):
        file_input=args.i[0]
    else:
        map_id=args.m[0]
    disease_id=args.d[0]

# Set the number of iterations
if args.s:
    repeats=args.s[0]
else:
    repeats=100000

read_line= []

#input file format: chr_name, coord_start, coord_end, quest_tag

cursor.execute("SELECT max(id) FROM enhancer_maps")
row=cursor.fetchone()[0]
if row is not None:
    map_id = row
    map_id= int(map_id)+1
else:
    map_id=1

name=basename(args.i[0]).split(".")[0]

#map_id=78
#name="T_Reg_Astrocytes_Brain_comb"

# Read new enhancer map and get the map id
map_id = load_map(file_input, name)     ### COMMENT THIS OUT FOR TESTING ###

# Count overlapping SNPs with enhancer regions from the new map for each leader and his follower sNPs 
main_pcount_hits(map_id, name)

# Calculate score for each disease based on the previously obtained hits
main_pcalculate_score(map_id, name)

# Significance testing - change the comments
main_calculate_significance(map_id, name, repeats, disease_id)
