#!/usr/bin/python 

#####################################################################################
# Script for counting the total number of leader and their follower SNPs associated #
# with  each disease in the DB. Called from integrateMyGWAS.py script.                #    
#####################################################################################

import MySQLdb as mdb
import sys
import csv
import re
from mySQLConnect import mysql_connect

con=mysql_connect()
cursor=con.cursor()

def count_SNPs():

    i=0

    curr_leader_name="0"
    curr_leader_name="0"
    curr_follower_name="0"

    num_leaders=0
    num_partial_leaders=0
    num_followers=0

    # Get all diseases from table DISEASE
    cursor.execute("SELECT id, name FROM diseases ORDER BY id");
    diseases =  cursor.fetchall()

    for disease in diseases:
        i+=1
    
        curr_disease_id = disease[0]
        curr_disease_name = disease[1]    

        # Count leaders
        cursor.execute("SELECT count(*) FROM leaders l, snps s WHERE s.id=l.snp_id AND disease REGEXP '^%s$|^%s,|,%s,|,%s$' AND EXISTS (SELECT * FROM leader_diseases WHERE snp_id=s.id AND disease_id=%s)", (curr_disease_id, curr_disease_id, curr_disease_id, curr_disease_id, curr_disease_id))
        num_leaders=cursor.fetchone()[0]
    
        # Count partial leaders 
        cursor.execute("SELECT count(*) FROM leaders l, snps s WHERE s.id=l.snp_id AND disease REGEXP '^%s$|^%s,|,%s,|,%s$' AND EXISTS (SELECT * FROM follower_diseases WHERE snp_id=s.id AND disease_id=%s)", (curr_disease_id, curr_disease_id, curr_disease_id, curr_disease_id, curr_disease_id))
        num_partial_leaders=cursor.fetchone()[0]

        # Count followers
        cursor.execute("SELECT count(*) FROM followers f, snps s WHERE s.id=f.snp_id AND disease REGEXP '^%s$|^%s,|,%s,|,%s$' AND EXISTS (SELECT * FROM follower_diseases WHERE snp_id=s.id AND disease_id=%s)", (curr_disease_id, curr_disease_id, curr_disease_id, curr_disease_id, curr_disease_id))
        num_followers=cursor.fetchone()[0]
    
        cursor.execute("UPDATE diseases SET num_leaders=%s, num_followers=%s WHERE id=%s", (num_leaders, num_followers+num_partial_leaders, curr_disease_id))
        con.commit()

    print "\nCounting done."
