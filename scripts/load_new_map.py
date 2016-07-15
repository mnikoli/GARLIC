#!/usr/bin/python

##############################################################################################################
# DELETES ALL ENHANCER_REGIONS FROM TABLE ENHANCER_REGIONS BEFORE SAVING A NEW REGULATORY MAP             
# A script for putting enhancer maps (format: chr start end) into tables enhancer_maps and enhancer_regions
# Parameters that need to be provided: path to a txt/bed file that contains enhancer intervals and a name 
# that is going to be used for that enhancer map.                                                            
##############################################################################################################

from mySQLConnect import mysql_connect
import MySQLdb as mdb
import sys
import csv
import re
import os
from os.path import basename
from remove_maps import removeMap

def load_map(file_input, name):
    con=mysql_connect()
    cursor=con.cursor()

    f=open(file_input)

    read_line= []
    
    try:
        cursor.execute("SELECT max(id) FROM enhancer_maps")
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)

    row=cursor.fetchone()[0]
    if row is not None:
        map_id = row
        map_id= int(map_id)+1
    else:
        map_id= 1
    
    # Get the parent dir
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


    # Create reports dir 
    results_dir = os.path.join(parent_dir, "reports")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Create testedCellTypes dir 
    results_dir = os.path.join(os.path.join(parent_dir, "reports"), "testedCellTypes")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # Create testedDiseases dir 
    results_dir = os.path.join(os.path.join(parent_dir, "reports"), "testedDiseases")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    try:
        cursor.execute("INSERT INTO enhancer_maps (name, id) VALUES(%s,%s)", (name, map_id))
        cursor.execute("DELETE FROM enhancer_regions");
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)

    curr_start=0
    curr_end=0
    k=0
    i=0
    for line in f:
        #...read line by line from the txt file...
        read_line.append(line.strip(''))
        row_el = re.split('\s+', read_line[k]) # re.split() used over str.split() because some lines have separator '\t' while some have one space ' '
        k=k+1
        
        #parse BED file so that chr_id corresponds to database chr_id
        chr_name=row_el[0]
        rex = re.match(r"chr(\w+)", chr_name)
        if rex is None:
            print "Error (Wrong format): ", line
            print "Input file should in BED format (at least TAB-delimited containing columns: chr_name, coord_start, coord_end)" 
            print "Exiting.."
            sys.exit(1)

        else:
            chr_id=rex.group(1)
        if chr_id:
            if chr_id == 'x' or chr_id =='X':
                chr_id=23
            if chr_id == 'y' or chr_id == 'Y':
                chr_id=24
            chr_id = str(chr_id)
            chr_id.strip()
            chr_id=int(chr_id)
    
        if int(row_el[1])<int(row_el[2]):
            curr_start=row_el[1]
            curr_start.strip()
            curr_start = int(curr_start)
            curr_end = row_el[2]
            curr_end.strip()
            curr_end = int(curr_end)
        else:
            print "Error (Wrong format): ", line
            print "Input file should be TAB-delimited containing at least columns: chr_name, coord_start, coord_end"    
            print "Exiting.."
            sys.exit(1)
    
        try:
            cursor.execute("INSERT INTO enhancer_regions (id, chr_id, start, end) VALUES(%s,%s,%s,%s)", (map_id, chr_id, curr_start, curr_end))
            cursor.execute("INSERT INTO stored_enhancer_regions (id, chr_id, start, end) VALUES(%s,%s,%s,%s)", (map_id, chr_id, curr_start, curr_end))
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)
    
    con.commit()

    print "\nMap id: ", map_id, "\nMap name: ", name, "\nTotal # of intervals read: ", k, "\n"
    f.close
    if k == 0:    
        print "Input file does not contain any regions, aborting."
        removeMap(map_id)
        sys.exit(1)
    return map_id
