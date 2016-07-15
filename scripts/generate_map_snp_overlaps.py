#!/usr/bin/python
###################################################################
# Script is used to dentify (leader and follower) SNPs associated 
# with a given disease that overlap a given regulatory  map.
###################################################################
    
import MySQLdb as mdb
import sys
import csv
import re
from os.path import basename
from mySQLConnect import mysql_connect
import argparse
from load_new_map import load_map
from remove_maps import removeMap

def checkChrId(x):
    chr_id=''
    if x==23:
        chr_id='X'
    elif x==24:
        chr_id='Y'
    else:
        chr_id=str(x)
    return chr_id

def mapSNPOverlaps(*args):
    disease_ids=[]
    remove_flag=0

    if ((args[0].m is None and args[0].i is None) or (args[0].m is not None and args[0].i is not None)):
        print "\nUsage (regulatory map already in database): garlic generateMapSNPOverlaps -m REG_MAP_ID <disease_id(s)>"
        print "Usage (regulatory map not in database): garlic generateMapSNPOverlaps -i REG_MAP_PATH <disease_id(s)> [-n REG_MAP_NAME]\n"
        print "If more than one disease ids are used with -d parameter, make sure to separate them with comma.\n" 
        print "Error: To few arguments"
        sys.exit(1)
    else:
        if args[0].i is not None:
            if args[0].n is None:
                base=basename(args[0].i[0]).split(".")
                name=base[0]
            else:
                name=args[0].n
    
            print "Saving new regulatory map in database..."
            curr_map_id=load_map(args[0].i[0], name)
            remove_flag=1
        else:
            curr_map_id=args[0].m[0]
    
        if args[0].d is not None:
            diseases=args[0].d[0]
            disease_list=diseases.split(',')
            for d in disease_list:
                disease_ids.append(int(d.strip()))
    
    con=mysql_connect()
    cursor=con.cursor()
    
    print "# Input parameters:"
    print "# Regulatory_map "+str(curr_map_id)
    if args[0].d is None:
        print "# Disease_ids: All diseases"
    else:
        print "# Disease_ids: "+str(disease_ids)
    print "#-------------------"
    
    disease_names = {}
    cursor.execute("SELECT id, name FROM diseases ORDER BY id")
    diseases = cursor.fetchall()
    for did in diseases:
        disease_names[int(did[0])]=str(did[1])
        if args[0].d is None:
            disease_ids.append(int(did[0]))
    read_line= []
    row_el= []
    
    for did in disease_ids:
        print "\nDisease: "+str(disease_names[did]), did
        database_snps = cursor.execute("SELECT s.name, start, chr_pos, end, l.chr_id, disease FROM leaders l, stored_enhancer_regions e, snps s WHERE s.id=l.snp_id AND (chr_pos BETWEEN start AND end) AND l.chr_id=e.chr_id AND e.id=%s AND disease REGEXP '^%s$|^%s,|,%s,|,%s$' ORDER BY (s.name)",(curr_map_id, did, did, did, did));
        
        i=0
        j=0
        if database_snps > 0:
            database_snps = cursor.fetchall()
            print "\nLeader SNPs:\n"    
            print "\tName\tOverlapping_region\tChr\tSNP_type\tAssoc_disease(s)"
    
            for snp in database_snps:
                split_disease=snp[5].split(",")
                disease_assoc_list=[]
                disease_assoc=''
                
                for sd in split_disease:
                    disease_assoc_list.append(disease_names[int(sd)])
                    #print disease_assoc_list
                    disease_assoc=str( ';'.join(str(x) for x in disease_assoc_list))

                print "\t"+str(snp[0])+"\t"+str(snp[1])+" < "+str(snp[2])+" < "+str(snp[3])+"\t"+"chr"+checkChrId(snp[4])+"\tL\t"+disease_assoc
                i+=1    
        
        database_snps = cursor.execute("SELECT s.name, start, chr_pos, end, f.chr_id, disease, following FROM followers f, stored_enhancer_regions e, snps s WHERE s.id=f.snp_id AND (chr_pos BETWEEN start AND end) AND f.chr_id=e.chr_id AND e.id=%s AND disease REGEXP '^%s$|^%s,|,%s,|,%s$' ORDER BY (s.name)",(curr_map_id, did, did, did, did));
        j=0
        if database_snps > 0:
            database_snps = cursor.fetchall()
            print "\n\nFollower SNPs:\n" 
            print "Followed_leader_SNP(s)\tName\tOverlapping_region\tChr\tSNP_type\tAssoc_disease(s)"
    
            for snp in database_snps:
                split_disease=snp[5].split(",")
                disease_assoc_list=[]
                disease_assoc=''
                
                split_following=snp[6].split(",")
                following_list=[]
                following=''
                
                for sf in split_following:
                    names = cursor.execute("SELECT name from snps WHERE id=%s" % sf);
                    if names>0:
                        name = cursor.fetchone()[0]
                        following_list.append(str(name.strip()))
                following = str(';'.join(str(x) for x in following_list))
                for sd in split_disease:
                    disease_assoc_list.append(disease_names[int(sd)])
                    disease_assoc=str( ';'.join(str(x) for x in disease_assoc_list))

                    print following+"\t"+str(snp[0])+"\t"+str(snp[1])+" < "+str(snp[2])+" < "+str(snp[3])+"\t"+"chr"+checkChrId(snp[4])+"\tF\t"+disease_assoc
                    j+=1
            
        print "\nOverlapping leader SNPs: "+str(i)
        print "Overlapping follower SNPs: "+str(j)
        print "In total: "+str(i+j)
        print "#---------------------------"
        
        if remove_flag:
            removeMap(curr_map_id)
