#!/usr/bin/python

#########################################################################################
# A script that parses the input from CSV file with GWAS data (published or unpublished) 
# using a wrap-up script integrateMyGWAS.py script                     
#
# Please refer to Manual regarding input format and available options. 
#########################################################################################

import MySQLdb as mdb
from mySQLConnect import mysql_connect
import sys
import csv
import re

con=mysql_connect()
cursor=con.cursor()
    
def load_gwas(gwas_input):

    k=0
    read_line= []
    row_el= []
    snp_risk=''

    f=open(gwas_input,'r')
    f.readline()

    pubmed_id=''
    disease_name=''
    chr_id=0
    chr_pos=0;
    mapped_gene=''
    strong_snp_risk=''
    snp_name=''
    p_val=''
    disease_gwas_ref=''
    k=0
    i=0
    print "Parsing GWAS input file "+str(gwas_input)

    for line in f:
        ### Read line by line from the csv file
        read_line.append(line.strip())
        row_el = read_line[k].split("\t")
        k=k+1
        pubmed_id='0'    
        if row_el[5] != 'NR' and row_el[5] != 'NS':
            rex = re.match(r"(rs[0-9]+)[-]*([-ACTG?,]*)",row_el[5])
            rex1 = re.search(r"([0-9]*\s*-\s*SNP|[Mm]arker\s*[Hh]aplotype\s*-\s*)",row_el[5])
            rex_snp_cluster = re.search(r",", row_el[6])

            pubmed_id=str(row_el[0])
            disease_name=row_el[1]
            if row_el[2] is '':
                chr_id=0
            else:
                chr_id=row_el[2]
            if row_el[3] is '':
                          chr_pos=0
            else:
                chr_pos=row_el[3]
                mapped_gene=row_el[4]
                p_val=row_el[7]

            ### Check whether current disease name already exists in the database and if not, save it in GARLIC DB
            if rex or rex_snp_cluster or rex1:
                cursor.execute("SELECT COUNT(*) FROM diseases WHERE name=%s" % disease_name);
                if cursor.fetchone()[0] == 0:    
                    print "Creating new entry for disease ", disease_name, "in table diseases"
                    cursor.execute("SELECT max(id) FROM diseases");
                    new_disease_id = cursor.fetchone()[0] 
                    if new_disease_id == None:
                        new_disease_id=1
                    else:
                        new_disease_id +=1
    
                    cursor.execute("INSERT INTO diseases (name, id) VALUES(%s, %s)", (disease_name, new_disease_id));
                    con.commit()
                    cursor.execute("SELECT id from diseases WHERE name=%s" % disease_name);
                    disease_id = cursor.fetchone()[0]
            
                snp_cluster = []
                snp_cluster = row_el[6].split(",")
    
                t=0
                for snp in snp_cluster:
                    
                    ### First, take only elements that have proper SNP ids
                    if re.match(r'\s*rs[0-9]+',snp):
                        t+=1
                        i+=1
                        snp_name = snp.strip()
                    
                        if t==1 and not rex1:
                            strong_snp_risk = rex.group(2)                
                        else:
                            strong_snp_risk = '0'
    
                        ### Next, check whether current SNP id already exists in the database and if not, create new rows in table "leaders" and "snps"
                        cursor.execute("SELECT COUNT(*) FROM snps WHERE name=%s" % snp_name);
                        if cursor.fetchone()[0] == 0:
                            print "Creating new entry in table leaders for snp ", snp_name
                
                            con.begin();
                            cursor.execute("SELECT max(id) FROM snps");
                            new_snp_id = cursor.fetchone()[0]
                            if new_snp_id == None:
                                new_snp_id = 1
                            else:
                                new_snp_id += 1
            
                            l=str(1)

                            cursor.execute("INSERT INTO leaders (snp_id, disease, chr_id, chr_pos, mapped_gene, strong_snp_risk, p_val, pubmed_id, disease_gwas_ref) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s)",(new_snp_id, disease_id, chr_id, chr_pos, mapped_gene, strong_snp_risk, p_val, pubmed_id, disease_id));
                            #print new_snp_id, snp_name, disease_id, chr_id, chr_pos, mapped_gene, strong_snp_risk, p_val, pubmed_id, disease_id
                            cursor.execute("INSERT INTO snps (id, name, leader) VALUES(%s,%s,%s)",(new_snp_id, snp_name, l));
                            cursor.execute("INSERT INTO leader_diseases (snp_id, chr_id, disease_id) VALUES(%s, %s, %s)",(new_snp_id, chr_id, disease_id)); 
                            con.commit()
                        else:
                            ### Leader SNP already exists, updationg DB info

                            # Get the id for current leader name 
                            cursor.execute("SELECT id FROM snps where name=%s" % snp_name);
                            curr_snp_id = cursor.fetchone()[0]
                        
                            # Get the list of previous snp-asociated diseases and strong_snp_risk nucleotides
                            print "Leader already exists, updating info for SNP "+str(snp_name)
                            cursor.execute("SELECT disease, strong_snp_risk, disease_gwas_ref FROM leaders WHERE snp_id=%s" % curr_snp_id);
                            row = cursor.fetchall()[0]
                            disease_list = row[0].split(",")
                            disease_list.append(str(disease_id))
                            new_disease_list = set(disease_list)
                            new_disease_list = str( ','.join(str(x) for x in new_disease_list))
                            new_snp_risk = row[1]
                            new_snp_risk = str(new_snp_risk).strip() + ","+ str(strong_snp_risk).strip()

                            # Reference list for gwas studies (how many gwas studies for this particular leader-SNP exist)
                            new_disease_gwas_ref = row[2]
                            new_disease_gwas_ref = str(new_disease_gwas_ref) + "," + str(disease_id)

                            # Update table leaders
                            cursor.execute("UPDATE leaders SET disease=%s, strong_snp_risk=%s, disease_gwas_ref=%s WHERE snp_id=%s", (new_disease_list, new_snp_risk, new_disease_gwas_ref, curr_snp_id));
                        
                            # Lastly, update table leader_diseases
                            cursor.execute("SELECT count(*) from leader_diseases where snp_id=%s and disease_id=%s", (curr_snp_id, disease_id));
                            check_id = cursor.fetchone()[0]
                            if check_id == 0:
                                cursor.execute("INSERT INTO leader_diseases (snp_id, chr_id, disease_id) VALUES(%s, %s, %s)",(curr_snp_id, chr_id, disease_id));

                            con.commit()
            else:
                print "Wrong input format", line
                sys.exit(1)
    
    print "\n\nSNPs read:", i
    print "Lines read:", k
    f.close
    print "Finished parsing GWAS input file\n\n"
