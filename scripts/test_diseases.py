#!/usr/bin/python
################################################################
# Main wrapper-script that calls all other scripts             #
# when testing diseases and traits with a given regulatory map #
# input file format: BED (chr, start, end) file (without chrM) # 
################################################################

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
from pcalculate_significance6 import main_pcalculate_significance
from calculate_dt_fdr import calculateDTFDR

con=mysql_connect()
cursor=con.cursor()

def testDiseases(*args): 
    if args[0].path is None:
        print "Error: Please provide a path to regulatory map.\n"
        print "Use -h option to get more information."
        sys.exit(1)
    else:
        file_input=args[0].path
        ltrshld=args[0].l

    if args[0].n is None or args[0].n == '':
        base=basename(args[0].path).split(".")
        name=base[0]
    else:
        name=args[0].n[0]

    repeats=args[0].s
    if args[0].t is None:
        ncpus=multiprocessing.cpu_count()
        print "Using "+str(multiprocessing.cpu_count())+" processes.\n"
    else: 
        ncpus=args[0].t
        if ncpus > multiprocessing.cpu_count():
            print "WARNING: You don't seem to have enough processes available on your machine."
            print "Using "+str(multiprocessing.cpu_count())+" processes instead of "+str(ncpus)+".\n"

    # Read new enhancer map and get the map id
    print "# Saving regulatory map in database..."
    map_id = load_map(file_input, name)
    
    print "# Analysing new regulatory map..."
    # Count overlapping SNPs with enhancer regions from the new map for each leader and his follower sNPs 
    main_pcount_hits(map_id, name)
    
    # Calculate score for each disease based on the previously obtained hits
    main_pcalculate_score(map_id, name)
    
    print "# Calculating disease/trait scores and empirical p values..." 
    # Calculate empirical p values
    main_pcalculate_significance(map_id, name, ltrshld, ncpus, repeats)

    # Recalculate FDR values for all diseases/traits
    calculateDTFDR()

    print "\n# Done."
