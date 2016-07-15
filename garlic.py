#!/usr/bin/python

####################################################
# Main wrapper-script that calls all other scripts 
####################################################
import sys, os

### Get the parent dir
parent_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(parent_dir, "scripts"))

import MySQLdb as mdb
from mySQLConnect import mysql_connect

from test_diseases import testDiseases
from test_cell_types import testCellTypes
from generate_region_snp_overlaps import regionSNPOverlaps
from generate_map_snp_overlaps import mapSNPOverlaps
from analyse_combinations import analyseCombinations
from integrate_gwas import integrateGWAS

import csv
import re
import argparse
from view_table_data import viewTableData
from os.path import basename
from load_new_map import load_map

from pcount_hits2 import main_pcount_hits
from pcalculate_score4 import main_pcalculate_score
from pcalculate_significance6 import main_pcalculate_significance

con=mysql_connect()
cursor=con.cursor()

# Parse the command line arguments
parser = argparse.ArgumentParser()
subparser = parser.add_subparsers(title='available parameters')

### Input parameters for different functions

# viewTableData
parser_viewTableData=subparser.add_parser('viewTableData')
parser_viewTableData.add_argument('-n', nargs=1, type=str, help="Database table name", metavar='TABLE_NAME')
parser_viewTableData.set_defaults(func=viewTableData)

# testDiseases
parser_testDiseases = subparser.add_parser('testDiseases')
parser_testDiseases.add_argument('path', type=str, help="Path to cell type or tissue regulatory map", metavar='path')
parser_testDiseases.add_argument('-n', default='', nargs=1, type=str, help="Given name for the new regulatory map", metavar='REG_MAP_NAME')
parser_testDiseases.add_argument('-l', default=7, nargs=1, type=int, help="Minimum number of GRRs in disease", metavar='MIN_GRR_NUM')
parser_testDiseases.add_argument('-t', default=8, nargs=1, type=int, help="Number of threads", metavar='NUM_THREADS')
parser_testDiseases.add_argument('-s', default=100000, nargs=1, type=int, help="Number of iterations for simulated score", metavar='NUM_ITERS')
parser_testDiseases.set_defaults(func=testDiseases)

# testCellTypes
parser_testCellTypes = subparser.add_parser('testCellTypes')
parser_testCellTypes.add_argument('did', type=int, help="Disease/trait id from database", metavar='did')
parser_testCellTypes.set_defaults(func=testCellTypes)

# integrateMyGWAS    
parser_regionSNPOverlaps = subparser.add_parser('integrateGWAS')
parser_regionSNPOverlaps.add_argument('gwas_input', type=str, help="Path to GWAS file", metavar='gwas_input')
parser_regionSNPOverlaps.add_argument('ld_input', type=str, help="Path to linkage-disequilibrium (LD) file", metavar='ld_input')
parser_regionSNPOverlaps.set_defaults(func=integrateGWAS)

# generateRegionSNPOverlaps
parser_regionSNPOverlaps = subparser.add_parser('generateRegionSNPOverlaps')
parser_regionSNPOverlaps.add_argument('c', type=str, help="Chromosome name", metavar='chr')
parser_regionSNPOverlaps.add_argument('s', type=int, help="Beginning of the region", metavar='start')
parser_regionSNPOverlaps.add_argument('e', type=int, help="End of the region", metavar='end')
parser_regionSNPOverlaps.set_defaults(func=regionSNPOverlaps)

# generateMapSNPOverlaps 
parser_mapSNPOverlaps = subparser.add_parser('generateMapSNPOverlaps')
parser_mapSNPOverlaps.add_argument('-d', nargs=1, type=str, help="One or more disease ids, separated by comma", metavar='did') 
parser_mapSNPOverlaps.add_argument('-m', nargs=1, type=str, help="Regulatory map id (if already in database)", metavar='REG_MAP_ID') 
parser_mapSNPOverlaps.add_argument('-i', nargs=1, type=str, help="Path to a regulatory map file in BED format (if not in database)", metavar='REG_MAP_PATH') 
parser_mapSNPOverlaps.add_argument('-n', nargs=1, type=str, help="Given name for the new regulatory map", metavar='REG_MAP_NAME') 
parser_mapSNPOverlaps.set_defaults(func=mapSNPOverlaps)

# analyseCombinations
parser_analyseComb = subparser.add_parser('analyseCombinations')
parser_analyseComb.add_argument('d', type=int, help="Disease/trait id from database", metavar='did')
parser_analyseComb.add_argument('-c', nargs=1, type=int, help="Maximum number of regulatory maps in combination. Default is 2.", metavar='NUM_COMB_EL')
parser_analyseComb.add_argument('-t', nargs=1, type=int, help="Number of threads.", metavar='NUM_THREADS')
parser_analyseComb.add_argument('-s', nargs=1, type=float, help="Parameter to limit the number of 'seed' regulatory maps which are use d as a basis to make combinations", metavar='SEED')
parser_analyseComb.add_argument('-n', nargs=1, type=int, help="Number of candidate combinations to test. Default is 3.", metavar='NUM_COMB_TEST')
parser_analyseComb.add_argument('-p', nargs=1, type=float, help="Use only those regulatory maps for which input disease has p value greater than a given parameter. Default is 0.001.", metavar='P_VALUE')
parser_analyseComb.add_argument('-m', nargs=1, type=int, help="Magnitude of p value improvement. Default is 5.", metavar='P_VALUE_IMPROVEMENT')
parser_analyseComb.set_defaults(func=analyseCombinations)
    
# Parse input arguments
args=parser.parse_args()
args.func(args)
