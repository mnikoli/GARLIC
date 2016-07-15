#!/usr/bin/python

#########################################################################
# Should be used before filtering, right after reading snps from gwas    
# haploreg repositories. Called from integrateMyGWAS.py script.        
# Remove leaders that have chr_pos=0 from tables:              
# leaders, snps, leader_diseases, follower_diseases              
# scripts that have to be run again after this:                
# count_lf, count_hits, calculate_enhancer_cov, calculate_significance   
#########################################################################

import MySQLdb as mdb
import sys
import csv
import re
from mySQLConnect import mysql_connect
import MySQLdb as mdb

con=mysql_connect()
cursor=con.cursor()

def remove_incomplete_data():
    i=0
    cursor.execute("SELECT s.id, s.name FROM leaders l, snps s WHERE l.snp_id=s.id AND l.chr_pos=0")
    leaders=cursor.fetchall()

    curr_leader_name=''
    curr_leader_id=0
    i=0
    for l in leaders:
        i+=1    
        curr_leader_name = l[0]
        curr_leader_id = l[1]
        print curr_leader_name, curr_leader_id
    
        print "Removing from table leaders"
        cursor.execute("DELETE FROM leaders WHERE snp_id=%s", curr_leader_id);
    
        print "Removing from table snps"
        cursor.execute("DELETE FROM snps WHERE id=%s", curr_leader_id);

        print "Removing from table leader_diseases"
        cursor.execute("DELETE FROM leader_diseases WHERE snp_id=%s", curr_leader_id);

        print "Removing from table follower_diseases"
        cursor.execute("DELETE FROM follower_diseases WHERE snp_id=%s", curr_leader_id);

        con.commit()
        print "Total number of leaders with unknown chromosome position removed:", i

    ##################################################
    #!      !       !       !       !       !       !#
    # Before removing SNPs with x in their name,     #
    # check with which disease are they associated.  # 
    # If this disease has 0 leaders, remove it.      #
    ##################################################

    cursor.execute("SELECT name, id FROM snps WHERE name regexp 'x' OR name regexp ':'")
    leaders=cursor.fetchall()

    curr_leader_name=''
    curr_leader_id=0
    j=0
    for l in leaders:
        j+=1
        curr_leader_name = l[0]
        curr_leader_id = l[1]
        print curr_leader_name, curr_leader_id
    
        print "Removing from table leaders"
        cursor.execute("DELETE FROM leaders WHERE snp_id=%s", curr_leader_id);
    
        print "Removing from table snps"
        cursor.execute("DELETE FROM snps WHERE id=%s", curr_leader_id);
    
        print "Removing from table leader_diseases"
        cursor.execute("DELETE FROM leader_diseases WHERE snp_id=%s", curr_leader_id);
    
        print "Removing from table follower_diseases"
        cursor.execute("DELETE FROM follower_diseases WHERE snp_id=%s", curr_leader_id);
    
        print "Removing from table genomic_regions"
        cursor.execute("DELETE FROM genomic_regions WHERE snp_id=%s", curr_leader_id);
        
    con.commit()
    print "Total number of leaders with x or : in the name:", j
