#/!/usr/bin/python 
#############################################################
# - PARALLELIZED VERSION                  
# - Count the number of SNPs (leaders and followers) occuring 
# within a certain enhancer map, without caring for SNP ids. 
#############################################################

import MySQLdb as mdb
import mySQLConnect
from mySQLConnect import mysql_connect
import sys
import csv
import re
import os
import threading
import multiprocessing
from multiprocessing import Process

def worker(l, curr_map_id, sem, curr_disease_id):
    p = multiprocessing.Process(target=count_hits, args=(l, curr_map_id, sem, curr_disease_id))
    p.start()
    p.join()
    sem.release()

def count_hits(l, curr_map_id, sem, curr_disease_id):
    hits=0
    curr_leader_name=l[0]
    curr_leader_id=l[1]

    con=mysql_connect()
    cursor=con.cursor()

    ### Perform counting only if the lSNP & disease_id combination does not exist in table HITS
    try:
        cursor.execute("SELECT COUNT(*) FROM hits WHERE snp_id=%s AND disease_id=%s AND map_id=%s", (curr_leader_id, curr_disease_id, curr_map_id))
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)

    if cursor.fetchone()[0] == 0:

        ### Check if the current lSNP overlaps with current DNAseI map
        try:
            cursor.execute("SELECT count(*) FROM leader_diseases lf WHERE lf.snp_id=%s AND lf.disease_id=%s AND EXISTS (SELECT * FROM enhancer_regions e, leaders l WHERE e.id=%s AND l.chr_pos BETWEEN e.start AND e.end AND l.snp_id=lf.snp_id and l.chr_id=lf.chr_id and l.chr_id=e.chr_id)",(curr_leader_id, curr_disease_id, curr_map_id))
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)
        leader_overlaps=cursor.fetchone()[0]

        ### Get followers from current leader
        try:
            cursor.execute("SELECT f.snp_id, f.chr_id, f.chr_pos FROM followers f, follower_diseases fd WHERE fd.snp_id=f.snp_id AND fd.disease_id=%s AND (f.following LIKE %s OR f.following LIKE %s OR f.following LIKE %s OR f.following LIKE %s) AND EXISTS (SELECT * FROM follower_diseases WHERE snp_id=f.snp_id AND disease_id=%s)", (curr_disease_id, str(curr_leader_id)+",%", "%,"+str(curr_leader_id), "%,"+str(curr_leader_id)+",%", curr_leader_id, curr_disease_id))
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)
        followers=cursor.fetchall()

        ### Get partial leaders following current leader
        try:
            cursor.execute("SELECT l.snp_id, l.chr_id, l.chr_pos FROM leaders l, follower_diseases fd WHERE fd.snp_id=l.snp_id AND fd.disease_id=%s AND (l.following LIKE %s OR l.following LIKE %s OR l.following LIKE %s OR l.following LIKE %s) AND EXISTS (SELECT * FROM follower_diseases WHERE snp_id=l.snp_id AND disease_id=%s)",(curr_disease_id, str(curr_leader_id)+",%", "%,"+str(curr_leader_id), "%,"+str(curr_leader_id)+",%", curr_leader_id, curr_disease_id))
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)

        pfollowers=cursor.fetchall()

        ### Create a unique list of followers for current leader
        all_followers = followers + pfollowers

        ### Count SNP overlaps with current DNAseI map from unique list of followers
        for af in all_followers:
            try:
                cursor.execute("SELECT count(*) from enhancer_regions WHERE id=%s AND chr_id=%s AND start<=%s AND end>=%s", (curr_map_id, af[1], af[2], af[2] ))
            except mdb.Error, e:
                print "Error %d: %s" % (e.args[0], e.args[1])
                sys.exit(1)
            followers_overlap=cursor.fetchone()[0]

            hits+=followers_overlap

        hits+=leader_overlaps

        try:
            cursor.execute("INSERT INTO hits (disease_id, snp_id, map_id, hits) VALUES(%s,%s,%s,%s)", (curr_disease_id, curr_leader_id, curr_map_id, hits))
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)
        con.commit()

def main_pcount_hits(map_id, name):
    con=mysql_connect()
    cursor=con.cursor()
    curr_map_id=map_id

    ### Get parent dir
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    ### Check for results dir 
    results_dir = os.path.join(parent_dir, "reports")
    if not os.path.exists(results_dir):
        print "Error: Directory "+results_dir+" does not exist!"
        sys.exit(1)
    
    ### Get all dieases from table DISEASE
    try:
        cursor.execute("SELECT id, name FROM diseases ORDER BY id");
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)
    diseases =  cursor.fetchall()

    i=0
    curr_leader_name="0"
    curr_follower_name="0"
    ncpus = multiprocessing.cpu_count()   # Number of threads to create

    ### Check whether a given map id exists in DB
    try:
        cursor.execute("SELECT COUNT(*) from enhancer_maps where id=%s" % curr_map_id)
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)
    if cursor.fetchone()[0]==0:
        print "Provided enhancer map ID does not exist in database"
        sys.exit(1)

    num_leaders=0
    num_followers=0
    hits=0
    
    ### There might be SNPs (in table leaders) that represent a leader for some diseases and for some they are denoted as follower SNPs.  (so-called  Partial leader.)
    ### That's why we need to iterate by diseases first
    for disease in diseases:
        i+=1
        curr_disease_id = disease[0]
        curr_disease_name = disease[1]    
        
        ### Get all leader SNPs from the current disease (table leader_diseases)
        try:
            cursor.execute("SELECT s.name, s.id FROM snps s, leader_diseases lf WHERE lf.disease_id=%s AND lf.snp_id=s.id" % curr_disease_id)
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)
        leaders=cursor.fetchall()
        
        workers = list()
        sem = multiprocessing.Semaphore(ncpus)
        ### For each leader SNP function count_hits is called in parallel
        for l in leaders:
            workers.append(threading.Thread(target=worker, args=(l, curr_map_id, sem, curr_disease_id)))
            sem.acquire()
            workers[-1].start()

        for w in workers: w.join()
