#!/usr/bin/python 
##############################################################################
# Script that calculates FDR for pvalues for each tested regulatory map  
##############################################################################

import sys
import csv
import os
import xlwt
from mySQLConnect import mysql_connect
import MySQLdb as mdb

from rpy2.robjects.packages import importr
stats = importr('stats', robject_translations={'format_perc': '_format_perc'})
from rpy2.robjects.vectors import FloatVector

from decimal import Decimal
from scipy import stats as s
import numpy

con=mysql_connect()
cursor=con.cursor()

def calculateDTFDR():
    init_num_set=[]
    generated_leader_ids=[]

    rows=0
    cols=0

    print "\nRecalculating FDR values for every disease/trait..."

    ### Get disease/trait details
    cursor.execute("SELECT name, id, num_leaders, num_followers FROM diseases WHERE EXISTS (SELECT * FROM results WHERE disease_id=id) ORDER BY id")
    details=cursor.fetchall()

    for info in details:
        curr_disease_id=info[1]
        curr_disease_name=info[0]

        ### Get all trait details from table results
        cursor.execute("SELECT name, pval, disease_score, count, disease_id, map_id FROM results, enhancer_maps WHERE disease_id=%s AND map_id=id ORDER BY pval" % curr_disease_id)    
        rows=cursor.fetchall()
        count=len(rows)

        ### Save pvalues in a vector and correct them with BH
        p_list=[]
        for row in rows:
            p_list.append(row[1])

        p_adjust=stats.p_adjust(FloatVector(p_list), method='BH')

        p=0
        for row in rows:
            cursor.execute("UPDATE results SET corrected_cell_pval=%s WHERE map_id=%s AND disease_id=%s", (round(p_adjust[p],10), row[5], row[4] ))
            p+=1

    con.commit()
