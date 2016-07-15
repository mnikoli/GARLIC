#!/usr/bin/python
############################################################################################
# Create haplotype blocks for each leader by going through his followers and adjusting min #
# and max values. Called from integrateMyGWAS.py script.                        #
############################################################################################

from mySQLConnect import mysql_connect
import sys
import csv
import re

con=mysql_connect()
cursor=con.cursor()

def define_GRRs():
    print "(Re)defining GRRs, this will take a while..."
    cursor.execute("TRUNCATE genomic_regions")
    curr_leader_name=''
    curr_leader_id=0
    curr_start=0
    curr_end=0
    curr_size=0
    curr_follower_name=''
    curr_follower_id=0
    curr_follower_pos=0
    curr_chr_id=0

    i=0
    j=0
    k=0

    #Get all diseases
    cursor.execute("SELECT id, name FROM diseases ORDER BY id")
    diseases = cursor.fetchall()

    # Foreach disease ...
    for d in diseases:
        k+=1
        curr_disease_id = d[0]
        curr_disease_name = d[1]

        i=0
        # ...Get all leaders assoc. with the current disease ...
        cursor.execute("SELECT s.name, id, chr_pos, chr_id FROM leaders l, snps s WHERE l.snp_id=s.id AND EXISTS (SELECT * FROM leader_diseases WHERE snp_id=s.id AND disease_id=%s)" % curr_disease_id)
        leaders=cursor.fetchall()

        # ...Foreach leader ...
        for l in leaders:
            followers=0
            partial_leaders=0    
            i+=1
            curr_leader_name = l[0]
            curr_leader_id = l[1]

            # Setup initial block coords by using chr_pos from current lSNP
            curr_chr_id = l[3]
            curr_start = l[2]
            curr_end = l[2]

            j=0
            # Get fSNPs from the current lSNP assoc. with the current disease
            cursor.execute("SELECT s.name, s.id, chr_pos, f.chr_id FROM followers f, snps s WHERE f.snp_id=s.id AND f.following REGEXP '^%s$|^%s,|,%s,|,%s$' AND EXISTS (SELECT * FROM follower_diseases WHERE snp_id=s.id AND disease_id=%s AND chr_id=f.chr_id)", (curr_leader_id, curr_leader_id, curr_leader_id, curr_leader_id, curr_disease_id))
            followers = cursor.fetchall()
    
            # Get followers from table leaders (partial leaders)
            cursor.execute("SELECT s.name, s.id, chr_pos, l.chr_id FROM leaders l, snps s WHERE l.snp_id=s.id AND l.following REGEXP '^%s$|^%s,|,%s,|,%s$' AND EXISTS (SELECT * FROM follower_diseases WHERE snp_id=s.id AND disease_id=%s AND chr_id=l.chr_id)", (curr_leader_id, curr_leader_id, curr_leader_id, curr_leader_id, curr_disease_id))
            partial_leaders=cursor.fetchall()
    
            all_followers = followers+partial_leaders
            for f in all_followers:
                j+=1
                curr_follower_name = f[0]
                curr_follower_id = f[1]
                curr_follower_pos = f[2]
#                print curr_follower_name, curr_follower_id, curr_follower_pos    
                if curr_follower_pos<curr_start:
                    curr_start = curr_follower_pos
                elif curr_follower_pos > curr_end:
                    curr_end = curr_follower_pos
            curr_size = curr_end - curr_start+1
            print curr_start, curr_end, curr_size
            cursor.execute("INSERT INTO genomic_regions (snp_id, disease_id, chr_id, startb, endb, sizeb, bin, num_followers) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)", (curr_leader_id, curr_disease_id, curr_chr_id, curr_start, curr_end, curr_size, 0, j));
    con.commit()
    
    print "Number of processed leaders: ", i
    print "Number of processed diseases: ", k
    print "\nDefining GRRs done."
