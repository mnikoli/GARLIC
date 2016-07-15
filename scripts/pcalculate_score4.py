#!/usr/bin/python
#################################################################         
# - Calculate the score for a given enhancer map by determining
# the coverage of enhancer regions for genomic region           
# represented by one leader SNP    
# - Paralelized version using multiprocessing library                       
#################################################################

### Here is assumed that DNAseI maps are in proper format (start<end) - checked in load_new_map.py

import MySQLdb as mdb
from mySQLConnect import mysql_connect
import sys
import csv
import os
import re
import threading
import multiprocessing
from multiprocessing import Process
import math

def worker(l, curr_map_id, sem):
    p = multiprocessing.Process(target=calculate_score, args=(l, curr_map_id, sem))
    p.start()
    p.join()
    sem.release()

def calculate_score(l, curr_map_id, sem):
    con=mysql_connect()
    cursor=con.cursor()

    curr_leader_id=0
    curr_startb=0
    curr_endb=0
    curr_sizeb=0
    chr_id=0

    num_followers=0
    ecov=0

    leader_score=0
    fe_score=0
    le_score=0
    re_score=0
    te_score=0
    ec_score=0

    curr_leader_id = l[0]
    curr_startb = l[1]
    curr_endb = l[2]
    curr_chr_id=l[3]
    num_followers=l[4]
    curr_sizeb=l[5]
    curr_disease_id=l[6]

    j=0
    
    ### Case I: Enhancer region contains the current genomic region (enhancer coverage is 100%)
    try:
        cursor.execute("SELECT start, end FROM enhancer_regions WHERE id=%s AND chr_id=%s AND start<=%s AND end>=%s", (curr_map_id, curr_chr_id, curr_startb, curr_endb))
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)
    enhancers_covering = cursor.fetchall()
    if enhancers_covering:
        ec_score=curr_endb - curr_startb + 1
        try:    
            cursor.execute("INSERT INTO overlapping_regions (snp_id, chr_id, start, end, disease_id, map_id, size) VALUES (%s,%s,%s,%s,%s,%s,%s)", (curr_leader_id, curr_chr_id, curr_startb, curr_endb, curr_disease_id, curr_map_id, curr_endb-curr_startb+1))
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)

        con.commit()
    else:
         ### CaseII: Whole enhancer regions are contained within the current genomic region
        try:
            cursor.execute("SELECT start, end FROM enhancer_regions WHERE id=%s AND chr_id=%s AND start>=%s AND start<=%s AND end<=%s", (curr_map_id, curr_chr_id, curr_startb, curr_endb, curr_endb))
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)
        
        enhancers_within = cursor.fetchall()

        for ew in enhancers_within:
            fe_score+=ew[1]-ew[0]+1

            try:
                cursor.execute("INSERT INTO overlapping_regions (snp_id, chr_id, start, end, disease_id, map_id, size) VALUES (%s,%s,%s,%s,%s,%s,%s)", (curr_leader_id, curr_chr_id, ew[0], ew[1]+1, curr_disease_id, curr_map_id, ew[1]-ew[0]+1))
            except mdb.Error, e:
                print "Error %d: %s" % (e.args[0], e.args[1])
                sys.exit(1)
        con.commit()
    
        ### Find all enhancer region(s) (if exist) that contain(s) left border of the haplotype block
        try:
            cursor.execute("SELECT start, end FROM enhancer_regions WHERE id=%s AND chr_id=%s AND start<=%s AND end>=%s AND end<=%s", (curr_map_id, curr_chr_id, curr_startb,curr_startb, curr_endb))
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)
        enhancers_lhb = cursor.fetchall()
        for elhb in enhancers_lhb:
            le_score+=elhb[1]-curr_startb+1
            try:                
                cursor.execute("INSERT INTO overlapping_regions (snp_id, chr_id, start, end, disease_id, map_id, size) VALUES (%s,%s,%s,%s,%s,%s,%s)", (curr_leader_id, curr_chr_id, curr_startb, elhb[1]+1, curr_disease_id, curr_map_id, elhb[1]-curr_startb+1))
            except mdb.Error, e:
                print "Error %d: %s" % (e.args[0], e.args[1])
                sys.exit(1)
 
        con.commit()

        ### Find all enhancer region(s) (if exist) that contain(s) right border of the haplotype block
        try:
            cursor.execute("SELECT start, end FROM enhancer_regions WHERE id=%s AND chr_id=%s AND start>=%s AND start<=%s AND end>=%s", (curr_map_id, curr_chr_id, curr_startb,curr_endb, curr_endb))
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)
        enhancers_rhb = cursor.fetchall()
        for erhb in enhancers_rhb:
            re_score+=curr_endb-erhb[0]+1
                        
            try:                    
                cursor.execute("INSERT INTO overlapping_regions (snp_id, chr_id, start, end, disease_id, map_id, size) VALUES (%s,%s,%s,%s,%s,%s,%s)", (curr_leader_id, curr_chr_id, erhb[0], curr_endb+1, curr_disease_id, curr_map_id, curr_endb-erhb[0]+1))
            except mdb.Error, e:
                print "Error %d: %s" % (e.args[0], e.args[1])
                sys.exit(1)
 
        con.commit()

    ecov=int(fe_score+le_score+re_score+ec_score)
    
    try:
        cursor.execute("UPDATE hits SET ecov=%s, score=%s, sizeb=%s WHERE snp_id=%s AND map_id=%s AND disease_id=%s", (ecov, ecov, curr_sizeb, curr_leader_id, curr_map_id, curr_disease_id))
        con.commit()
    except mdb.Error, e:
        if con:
            con.rollback()
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)

def main_pcalculate_score(map_id, name):
    con=mysql_connect()
    cursor=con.cursor()

    curr_map_id=map_id

    i=0
    s=0
    score=0
    fe_score=0 # full enhancer score
    le_score=0 # left haplotype border score within an enhancer 
    re_score=0 # right haplotype border score within an enhancer 

    # Check whether a given map id exists in DB
    try:
        cursor.execute("SELECT COUNT(*) from enhancer_maps where id=%s" % curr_map_id)
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)
    if cursor.fetchone()[0]==0:
        print "Provided enhancer map ID does not exist in database"
        sys.exit(1)
    
    # Get coordinates from genomic regions represented by one lSNP and one disease
    try:
        cursor.execute("SELECT distinct snp_id, startb, endb, chr_id, num_followers, sizeb, disease_id FROM genomic_regions")
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)
    leaders=cursor.fetchall()

    ncpus = multiprocessing.cpu_count()    #-1   # Number of threads to create
    
    workers = list()
    sem = multiprocessing.Semaphore(ncpus)

    ### For each genomic region, calculate the number of overlappng nts. in parallel     
    for l in leaders:
        workers.append(threading.Thread(target=worker, args=(l, curr_map_id, sem)))
        sem.acquire()
        workers[-1].start()
    for w in workers: w.join()
