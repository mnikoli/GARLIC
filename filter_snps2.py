#!/usr/bin/python

##########################################################################################################
# A script that identifies duplicate leaders (markers of common human haplotypes) and removes them in an #
# appropriate manner, since it has to deal with followers, SNPs that are in linkage-disequilibrium (LD)  #
# with duplicate leaders. Called from integrateMyGWAS.py script.                     #
##########################################################################################################

import MySQLdb as mdb
from mySQLConnect import mysql_connect
import sys
import csv
import re

con=mysql_connect()
cursor=con.cursor()

def filter_dhlSNPs():

    # Get all dieases from table DISEASE
    cursor.execute("SELECT id, name FROM diseases ORDER BY id");
    diseases =  cursor.fetchall()

    i=0
    j=0
    k=0
    n=0
    curr_disease_id = 0
    curr_disease_name = "0"
    curr_chr_id = 0
    curr_leader_name = "0"
    curr_leader_id = 0
    del_leader=0    # Flag - currently identified duplicate leader should be deleted from table leaders (for the current disease)
    del_follower=0
    count=0         # Flag that is used to mark the last leader that will represent a haplotype

    # For each disease..
    for disease in diseases:
        i+=1
    
        curr_disease_id = disease[0]
        curr_disease_name = disease[1]
        print "Disease id:"+str(curr_disease_id)+" "+curr_disease_name
    
        cursor.execute("SELECT DISTINCT chr_id FROM leaders WHERE disease LIKE %s OR disease LIKE %s OR disease LIKE %s OR disease LIKE %s GROUP BY(chr_id)", (str(curr_disease_id)+",%", "%,"+str(curr_disease_id), "%,"+str(curr_disease_id)+",%", curr_disease_id));
        chr_ids = cursor.fetchall()
    
        # ... And for each chromosome
        for chr_id in chr_ids:
            curr_chr_id = chr_id[0];
    
            ### Find leader SNPs in LD associated with the same disease (curr_disease_id) on the same chromosome (curr_chr_id) ###
            # Get all leaders (and their info) that are on the same chromosome & associated with the same disease
            cursor.execute("SELECT s.name, s.id, chr_pos, following, disease FROM leaders l, snps s WHERE l.snp_id=s.id AND s.leader=%s AND (disease LIKE %s OR disease LIKE %s OR disease LIKE %s OR disease LIKE %s) AND EXISTS(SELECT * FROM leaders WHERE snp_id=l.snp_id AND chr_id=l.chr_id AND disease=l.disease and chr_id=%s) ORDER by s.name", (1, str(curr_disease_id)+",%", "%,"+str(curr_disease_id), "%,"+str(curr_disease_id)+",%", curr_disease_id, curr_chr_id))

            leaders=cursor.fetchall()
            k=0
            m=0
        
            # For each leader...
            for l in leaders:
                k+=1
                curr_leader_name = l[0]
                curr_leader_id = l[1]
            
                llist_tmp=[]
                llist_tmp.append(str(curr_leader_id))
                                        
                # Check for duplicate haplotype leader(s) by looking at the following field to see whether they are in LD for the same disease
                cursor.execute("SELECT s.name, s.id, chr_pos, following, disease, r2, d, ref, alt, AFR, AMR, ASN, EUR, mapped_gene FROM leaders l, snps s WHERE l.snp_id=s.id AND s.leader=%s AND (disease LIKE %s OR disease LIKE %s OR disease LIKE %s OR disease LIKE %s) AND EXISTS(SELECT * FROM leaders WHERE snp_id=l.snp_id AND chr_id=l.chr_id AND disease=l.disease and chr_id=%s) ORDER BY s.name", (1, str(curr_disease_id)+",%", "%,"+str(curr_disease_id), "%,"+str(curr_disease_id)+",%", curr_disease_id, curr_chr_id))

                leaders1=cursor.fetchall()
                if len(leaders1)>=1:
                    m=0
                    for l1 in leaders1:
                        m+=1
                        if curr_leader_id!=l1[1] and str(l1[3])!='0': 
                            dhl_name = l1[0]
                            dhl_id = l1[1]
                            dhl_pos = l1[2]
                            dhl_following = str(l1[3])
                            dhl_following_list = dhl_following.split(",")
    
                            if set(dhl_following_list) & set(llist_tmp):
    
                                # ...Duplicate candidate leader identified... 
                                ### Check whether this leader is the real candidate or he is not a candidate anymore (data already filtered from relation to another leader that was in LD with this one) ###
                                # Get refreshed disease set & following set from duplicate candidate leader
                                cursor.execute("SELECT disease, following FROM leaders WHERE snp_id=%s", dhl_id)
                                refreshed_dhl_data = cursor.fetchall()[0]
                                dhl_diseases = refreshed_dhl_data[0]
                                new_dhl_following = refreshed_dhl_data[1]
                                new_dhl_following_int = new_dhl_following.split(",")
    
                                dhl_diseases_int = dhl_diseases.split(",")    
                                dhl_diseases_int2 = map(int, dhl_diseases_int)    
                                if set(new_dhl_following_int)&set(llist_tmp):    # No else: statement
                                    # It is a real duplicate leader for the same haplotype
                                    del_leader=0
                            
                                    # ...Check with how many diseases marked duplicate leader is associated  
                                    if len(dhl_diseases_int)>1:
                                        # Check whether other diseases have lower IDs than the current one... #
                                        # If true, data is already corrected for them and we can treat this   #        
                                        # leader as it has only one disease...                      #    

                                        still_a_leader=0
                                        if (max(dhl_diseases_int2) == int(curr_disease_id)):
                                            for dhld in dhl_diseases_int:
                                                cursor.execute("SELECT count(*) FROM leader_diseases WHERE snp_id=%s AND disease_id=%s", (dhl_id, dhld))
                                                disease_num=cursor.fetchone()[0]
                                                if disease_num > 0 and int(dhld) < max(dhl_diseases_int2):
                                                    still_a_leader=1

                                            if still_a_leader==0:    
                                                del_leader=1    
                                        
                                    elif len(dhl_diseases_int)==1:    
                                        del_leader=1
                                    else:
                                        print "Error in table leaders: disease field is empty!";
                                        sys.exit(1)

                                    # Ex-leader is now associated with current disease as a follower and not as a  #
                                    # leader anymore. Make changes in tables follower_diseases and leader_diseases #
                                    # Next, remove ex-leader's id from all other leaders that might be following him. 
                                    # Remove ex-leader's ID from other leaders that might be following him, only for the current disease! #
                                    if del_leader==1: 
                                        cursor.execute("SELECT s.name, l.following, s.id FROM leaders l, snps s WHERE (following LIKE %s OR following LIKE %s OR following LIKE %s OR following LIKE %s) AND s.id=l.snp_id", (str(dhl_id)+",%", "%,"+str(dhl_id), "%,"+str(dhl_id)+",%", dhl_id))
                                    else:    # Get all leaders that are following dhl (just for the current disease)
                                        cursor.execute("SELECT s.name, l.following, s.id FROM leaders l, snps s WHERE (following LIKE %s OR following LIKE %s OR following LIKE %s OR following LIKE %s) AND s.id=l.snp_id AND EXISTS (SELECT * FROM leader_diseases WHERE snp_id=s.id AND disease_id=%s and chr_id=l.chr_id)", (str(dhl_id)+",%", "%,"+str(dhl_id), "%,"+str(dhl_id)+",%", dhl_id, curr_disease_id))

                                    leaders_for_update=cursor.fetchall()
                                    print "Leaders for update: ", leaders_for_update, "dhl_id: ",dhl_id
                                    for lfu in leaders_for_update:
                                        lfu_name=lfu[0]
                                        lfu_following=lfu[1]
                                        lfu_id=lfu[2]
                                        print "Current lfu cmp: ", lfu_name, lfu_following, lfu_id
                                        if lfu_id!=dhl_id:
                                            lfu_following=lfu_following.split(",")
                                            tmp_following=[]
                                            print lfu_id, "!=",dhl_id
                                            # Remove id
                                            for u in lfu_following:
                                                
                                                tmp_following.append(int(u))
                                            
                                            tmp_following.remove(dhl_id)

                                            # If after the removal, no values are left - set value to 0
                                            if len(tmp_following)==0:
                                                tmp_following.append(0)
                                
                                            lfu_following=str( ','.join(str(x) for x in tmp_following))

                                            # Update column following in table leaders 
                                            cursor.execute("UPDATE leaders SET following=%s WHERE snp_id=%s", (lfu_following, lfu_id))
                                    con.commit() 

                                    ### Duplicate leader is marked for removal from table leaders into table followers) FOR THE CURRENT DISEASE ###
                                    # del_leader=1 - dhl is either assoc. with only 1 trait, or curr.disease has the max id in his disease field (last disease to be dealt with)
                                    if del_leader==1: 
                                        # Perform number of necessary changes in database            
                                        cursor.execute("SELECT COUNT(*) FROM followers WHERE snp_id=%s", (dhl_id))
                                        if cursor.fetchone()[0] == 0:
                                            # First, create the new disease set by using disease information from other leaders being followed by this ex-leader 
                                            disease_add=0
                                            new_disease_set=0
                                            dhl_diseases_new=dhl_diseases    

                                            # Get and add disease sets from other leaders that are being followed by ex-leader 
                                            for dlf in new_dhl_following_int:
                                                cursor.execute("SELECT s.name, disease FROM leaders l, snps s WHERE l.snp_id=s.id and s.id=%s",dlf)
                                                nds=cursor.fetchone()
                                                new_disease_set=nds[1]
                                                dhl_diseases_new=dhl_diseases_new+","+str(new_disease_set)

                                            dhl_diseases_new_int=dhl_diseases_new.split(",")
                                            dhl_diseases_new_int=set(dhl_diseases_new_int)
                                            dhl_diseases_new = str( ','.join(str(x) for x in dhl_diseases_new_int))
                                            print "Moving completed."    
                                             
                                            ### Then, move an entry from table leaders to table followers ###
                                            con.begin()
                                            # Get refreshed data from table leaders
                                            cursor.execute("SELECT r2, d, ref, alt, AFR, AMR, ASN, EUR, mapped_gene FROM leaders WHERE snp_id=%s", dhl_id)
                                            l2=cursor.fetchall()[0]

                                            # Create entry in table followers
                                            cursor.execute("INSERT INTO followers (snp_id, disease, chr_id, chr_pos, mapped_gene, following, r2, d, ref, alt, AFR, AMR, ASN, EUR) VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", (int(dhl_id), dhl_diseases_new, curr_chr_id, int(dhl_pos),l2[8], new_dhl_following, float(l2[0]), float(l2[1]), l2[2], l2[3], float(l2[4]), float(l2[5]), float(l2[6]), float(l2[7])))

                                            # Delete entry from table leaders
                                            cursor.execute("DELETE FROM leaders WHERE snp_id=%s", dhl_id)
    
                                            # Update entry in table snps - SNP is no longer a leader, but a follower FOR THE CURRENT DISEASE
                                            cursor.execute("UPDATE snps SET leader=%s WHERE id=%s", (0, dhl_id))
                                            con.commit()

                                    ### Next, move entry from table leader_diseases to follower_diseases for current disease ### (regardles whether del_leader is 1 or 0)
                                    cursor.execute("SELECT COUNT(*) FROM follower_diseases WHERE snp_id=%s AND disease_id=%s AND chr_id=%s", (dhl_id, curr_disease_id, curr_chr_id))
                                    if cursor.fetchone()[0] == 0:
                                        con.begin()
                                        # Insert entry into follower_diseases...
                                        cursor.execute("INSERT INTO follower_diseases (snp_id, chr_id, disease_id) VALUES(%s,%s,%s)", (dhl_id, curr_chr_id, curr_disease_id));
                                        # Delete same entry from leader_diseases
                                        cursor.execute("DELETE FROM leader_diseases WHERE snp_id=%s AND chr_id=%s AND disease_id=%s", (dhl_id, curr_chr_id, curr_disease_id))
                                        con.commit()

                                    ### Followers from the ex-leader need to be taken care of as well ###
                                    n=0
                                    cursor.execute("SELECT s.name, s.id, disease, following FROM followers f, snps s WHERE (following LIKE %s OR following LIKE %s OR following LIKE %s OR following LIKE %s) AND f.snp_id=s.id AND s.leader=0 AND EXISTS (SELECT * FROM follower_diseases WHERE s.id=snp_id and disease_id=%s) GROUP by s.name", (str(dhl_id)+",%", "%,"+str(dhl_id), "%,"+str(dhl_id)+",%", dhl_id, curr_disease_id))
                                    followers = cursor.fetchall()

                                    # For each follower of ex-leader check if
                                    for f in followers:
                                        n+=1
                                        curr_follower_name = f[0]
                                        curr_follower_id = f[1]
                                        curr_follower_disease_set = f[2]
                                        curr_follower_following = f[3]    
                                
                                        # he is following more than one leader
                                        curr_follower_following = curr_follower_following.split(",")
                                        if len(curr_follower_following)>1: 
                                            # Remove ex-leader's id from a follower with multiple (other) leaders 
                                            tmp_following=[]
                                            for c in curr_follower_following:
                                                tmp_following.append(int(c))
                                        
                                            tmp_following.remove(dhl_id)
                                            curr_follower_following= str( ','.join(str(x) for x in tmp_following)) # List of leaders for the current fSNP without the ex_lSNP
                                        
                                            curr_follower_disease_int = curr_follower_disease_set.split(",")
                                            intersect=list(set(curr_follower_disease_int)&set(dhl_diseases_int))
    
                                            # Remove current disease id from the intersected set
                                            intersect.remove(str(curr_disease_id))
                                            if intersect:
                                                # For each disease common to dhl and a current_follower (except the current disease)...
                                                for dh_dis in intersect:
                                                    # ... For each leader that the current SNP follower is following (except the ex-leader) ...
                                                    disease_found=0
                                                    for s in tmp_following:
                                                        # Look for a leader that is associated with the current disease (from the intersected set) and currently being followed by the current fSNP
                                                        r=cursor.execute("SELECT count(*) FROM leaders l, snps s WHERE l.snp_id=s.id AND s.id=%s AND (disease LIKE %s OR disease LIKE %s OR disease LIKE %s OR disease LIKE %s) AND exists (SELECT * FROM leader_diseases WHERE l.snp_id=snp_id and disease_id=%s)", (s, str(dh_dis)+",%", "%,"+str(dh_dis), "%,"+str(dh_dis)+",%", dh_dis, dh_dis))
                                                        if r>0:
                                                            counts=cursor.fetchone()[0]
                                                            if counts > 0:                
                                                                disease_found=1
                                                                break
                                                    if disease_found == 0 and del_leader == 1:
                                                        # Remove curr_disease_id from current followers' disease set
                                                        curr_follower_disease_int.remove(dh_dis)
                                                        new_curr_follower_disease_set = str( ','.join(str(x) for x in curr_follower_disease_int))

                                                        # Delete from table follower_diseases and update columns following and disease in table followers 
                                                        con.begin()
                                                        cursor.execute("DELETE FROM follower_diseases WHERE snp_id=%s AND disease_id=%s", (curr_follower_id, dh_dis))
                                                        cursor.execute("UPDATE followers SET following=%s, disease=%s WHERE snp_id=%s", (curr_follower_following, new_curr_follower_disease_set, curr_follower_id))
                                                        con.commit()
    
                                                    disease_found=0
                                            else:    
                                                cursor.execute("UPDATE followers SET following=%s WHERE snp_id=%s", (str( ','.join(str(x) for x in tmp_following)), curr_follower_id))
                                                con.commit()

                                        ### Follower SNP has one leader - ex-leader ###
                                        elif len(curr_follower_following)==1: 
                                            if del_leader==1: # If dhl is moved from table leaders to table followers                        
                                                # Update info for fSNP - set a new leader to be followered instead of the dhl, leave disease set intacked 
                                                cursor.execute("UPDATE followers SET following=%s WHERE snp_id=%s", (curr_leader_id, curr_follower_id))
                                            else:
                                                cursor.execute("UPDATE followers SET following=%s WHERE snp_id=%s", (str(curr_leader_id)+","+str( ','.join(str(x) for x in curr_follower_following)), curr_follower_id))
                                            con.commit()
                                        else:
                                            print "Error table followers: invalid field 'following'"
                                            print "fSNP: "+curr_follower_name
    print "\nSNP data ajdusting done."
