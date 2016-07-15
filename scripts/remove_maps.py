#!/usr/bin/python

#####################################################################
# Script that removes data from all tables for a set of DNAseI maps #    
#####################################################################

import MySQLdb as mdb
from mySQLConnect import mysql_connect
from calculate_dt_fdr import calculateDTFDR

import sys
import os

def removeMap(mid):
    re_fdr_flag=0
    con=mysql_connect()
    cursor=con.cursor()

    # Check if map mid exists in DB
    if cursor.execute("SELECT count(*) FROM enhancer_maps WHERE id=%s" % mid):    
        flag = cursor.fetchone()[0]
        if flag:
            # Check if fdr values should be re-calculated after regulatory map removal
            cursor.execute("SELECT count(*) FROM results WHERE map_id=%s" % mid)
            num_rows=cursor.fetchone()[0]
            if num_rows>1:
                re_fdr_flag=1
            # Get regulatory map name
            cursor.execute("SELECT name FROM enhancer_maps WHERE id=%s" % mid)
            name=cursor.fetchone()[0]
            print "Removing regulatory map %s from database..." % name

            try:
                cursor.execute("DELETE FROM enhancer_maps WHERE id=%s" % mid)
                cursor.execute("DELETE FROM stored_enhancer_regions WHERE id=%s" % mid)
                cursor.execute("DELETE FROM enhancer_regions WHERE id=%s" % mid)
                cursor.execute("DELETE FROM hits WHERE map_id=%s" % mid)
                cursor.execute("DELETE FROM results WHERE map_id=%s" % mid)
                cursor.execute("DELETE FROM overlapping_regions WHERE map_id=%s" % mid)
            except mdb.Error, e:
                print "Error %d: %s" % (e.args[0], e.args[1])
                sys.exit(1)

            con.commit()
            
            if re_fdr_flag:
                print "Recalculating FDR values"
                calculateDTFDR()

        else:
            print "Map with id %s does not exist in DB!" % mid
