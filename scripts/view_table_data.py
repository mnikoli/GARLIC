#!/usr/bin/python
#######################################################################
# A simple script that prints contents from table of interest from DB #
#######################################################################

import MySQLdb as mdb
import sys
import csv
import re
from mySQLConnect import mysql_connect
import argparse

con=mysql_connect()
cursor=con.cursor()

def viewTableData(*args):
    if args[0].n is None:
        cursor.execute("SHOW TABLES");
        tables = cursor.fetchall()
        print "Warning: Table name is missing, listing available tables:\n"
        for table in tables:
            print table[0]
        
        print "\nUse -h option to get more information."
        sys.exit(1)
    else:
        table_name=args[0].n[0]
        
        #Get column names
        try:
            cursor.execute("SHOW COLUMNS FROM %s" % table_name)
        except:
            print "Error: Provided table does not exist in database."
            print "\nUse -h option to get more information."
            sys.exit(1)

        print "# Table name: ", table_name
        print "#-------------------"

        cols=cursor.fetchall()
        colnames=""
        for cn in cols:
            colnames=colnames+str(cn[0])+"\t"
    
        print colnames.strip('\t')
        # Get rows
        cursor.execute("SELECT * FROM %s" % table_name)
        rows=cursor.fetchall()
        for row in rows:
            line=str( '\t'.join(str(r) for r in row))
            print line
