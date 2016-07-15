#!/usr/bin/python

########################################################
# Assign defined GRRs to bins according to their size. #
# Called from integrateMyGWAS.py script               #
########################################################

import MySQLdb as mdb
from mySQLConnect import mysql_connect
import sys
import os
import re

con=mysql_connect()
cursor=con.cursor()

def bin_GRRs():
    cursor.execute("SELECT MAX(sizeb) FROM genomic_regions")
    max_size=cursor.fetchone()[0]

    start_list=[1, 1001, 10001, 25001, 55000]
    end_list=[1000, 10000, 25000, 55000, max_size]

    for i in range(0, len(start_list)):
        print start_list[i], end_list[i]
        bin_id=i+1
        cursor.execute("UPDATE genomic_regions SET bin=%s WHERE sizeb BETWEEN %s AND %s", (bin_id, start_list[i], end_list[i]))
        con.commit()

    cursor.execute("UPDATE genomic_regions SET startb=startb-500+sizeb/2, endb=endb+500-sizeb/2, sizeb=endb-startb+1 WHERE sizeb<1000");
    con.commit()

    print "\nAssigning GRRs to bins done."
