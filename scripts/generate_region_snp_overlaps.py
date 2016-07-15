#!/usr/bin/python
#####################################################################
# Identify (leader and follower) SNPs that overlap the given locus. #
#####################################################################
    
import MySQLdb as mdb
import sys
import csv
import re
from mySQLConnect import mysql_connect
    
con=mysql_connect()
cursor=con.cursor()
    
def regionSNPOverlaps(*args): 
    print args[0].c, args[0].s, args[0].e
    chr_id_tmp=''
    disease_ids=[]
    if args[0].c is None or args[0].s is None or args[0].e is None:
        print "\nUsage: garlic generateRegionSNPOverlaps <chr> <start> <end>\n"
        print "Error: To few arguments"
        sys.exit(1)
    else:
        chr_id =args[0].c
        rex = re.match(r"chr(\w+)", chr_id)
        if rex is None:
            print "Error (Wrong format): ", chr_id
            sys.exit(1)
 
        else:
            chr_id=rex.group(1)
            if chr_id:
                if chr_id == 'x' or chr_id =='X':
                    chr_id_tmp='X'
                    chr_id=23
                elif chr_id == 'y' or chr_id == 'Y':
                    chr_id_tmp='Y'
                    chr_id=24
                elif int(chr_id)>=1 and int(chr_id)<=22:
                    chr_id = str(chr_id)
                    chr_id.strip()
                    chr_id=int(chr_id)
                    chr_id_tmp=chr_id

                else:
                    print "Error (Wrong format): chr"+str(chr_id)
                    sys.exit(1)

        start=args[0].s
        end=args[0].e
    
    print "\n# Query region: chr"+str(chr_id_tmp)+" "+str(start)+"-"+str(end)
    print "#-------------------"
    
    disease_names = {}
    cursor.execute("SELECT id, name FROM diseases ORDER BY id")
    diseases = cursor.fetchall()
    for did in diseases:
        disease_names[int(did[0])]=str(did[1])
    
    read_line= []
    row_el= []
    
#database_snps = cursor.execute("SELECT s.name, chr_pos, disease FROM leaders l, snps s WHERE s.id=l.snp_id AND l.chr_id=%s AND (chr_pos BETWEEN %s AND %s) GROUP BY (s.name)",(chr_id, start, end));
    database_snps = cursor.execute("SELECT s.name, chr_pos, disease FROM leaders l, snps s WHERE s.id=l.snp_id AND l.chr_id=%s AND (chr_pos BETWEEN %s AND %s) ORDER BY (s.name)",(chr_id, start, end));

        
    i=0
    if database_snps > 0:
        database_snps = cursor.fetchall()
        print "\nLeader SNPs:\n"    
        print "\tName\tChr\tOverlapping_region\tSNP_type\tAssoc_disease(s)"
    
        for snp in database_snps:
            split_disease=snp[2].split(",")
            disease_assoc_list=[]
            disease_assoc=''
            for sd in split_disease:
                disease_assoc_list.append(disease_names[int(sd)])
                #print disease_assoc_list
                disease_assoc=str( ';'.join(str(x) for x in disease_assoc_list))
            print "\t"+str(snp[0])+"\tchr"+str(chr_id_tmp)+"\t"+str(start)+" < "+str(snp[1])+" < "+str(end)+"\tL\t"+disease_assoc
            i+=1    
        
    database_snps = cursor.execute("SELECT s.name, chr_pos, disease, following FROM followers f, snps s WHERE s.id=f.snp_id AND f.chr_id=%s AND (chr_pos BETWEEN %s AND %s) ORDER BY s.name",(str(chr_id), start, end));
    j=0
    if database_snps > 0:
        database_snps = cursor.fetchall()
        print "\n\nFollower SNPs:\n" 
        print "Followed_leader_SNP(s)\tName\tChr\tOverlapping_region\tSNP_type\tAssoc_disease(s)"
    
        for snp in database_snps:
            split_disease=snp[2].split(",")
            disease_assoc_list=[]
            disease_assoc=''
                
            split_following=snp[3].split(",")
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
                    #print disease_assoc_list
                    disease_assoc=str( ';'.join(str(x) for x in disease_assoc_list))
            print following+"\t"+str(snp[0])+"\tchr"+str(chr_id_tmp)+"\t"+str(start)+" < "+str(snp[1])+" < "+str(end)+"\tF\t"+disease_assoc
            j+=1
    print "\nNumber of overlapping leader SNPs: "+str(i)
    print "Number of overlapping follower SNPs: "+str(j)
    print "Total number of overlapping SNPs: "+str(i+j)
    print "#---------------------------"
