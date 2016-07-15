#!/usr/bin/python 

####################################################################################
#               ***Parallelized as a pool of workers ***              
# A script that calculates "in silico" which combination of enhancer maps could    
# potentially increase the significance of a particular / all disease/s           
#                                                                                  
# This script prioritizes combinations by fold change first and then ranks them    
# remaining maps by their fold change                                               
####################################################################################

import MySQLdb as mdb
import sys
import csv
import re
import operator
from mySQLConnect import mysql_connect
import itertools

import threading
import multiprocessing
from multiprocessing import Process, Value
from functools import partial
import argparse

from subprocess import Popen, PIPE, call
import os
import ntpath
from remove_maps import removeMap
import subprocess

con=mysql_connect()
cursor=con.cursor()

def analyseOneComb(combinations_dict, covered_grs_dict, curr_scores_dict, input_disease_id, num_leaders,zero_ecov_list, non_existant, f):
    con=mysql_connect()
    cursor=con.cursor()
    
    tuple_dict={} # for each of map combinations, keys are chr ids 
    
    ### Initiate dict of tuples
    for t in range(1,25):
        tuple_dict[t]=[]

    tmp_list=combinations_dict[f]
    #print tmp_list
    snp_list_in_comb=[]
    id_list=[]
    
    ### Get overlapping details for each map from the current combination
    for el in tmp_list:
        cursor.execute("SELECT distinct chr_id, start, end, snp_id FROM overlapping_regions WHERE map_id=%s AND disease_id=%s", (el, input_disease_id))
        tuples=cursor.fetchall()
        for t in tuples:
            tuple_dict[t[0]].append((t[1], t[2]))
            snp_list_in_comb.append(t[3])

        id_list.append(curr_scores_dict[el][0])

    uniq_tuple_dict={}
    uniq_tuple_dict=merge_intervals(tuple_dict)
    score=sum_nts(uniq_tuple_dict)

    covered=num_leaders[1]-len(list(set(zero_ecov_list) - (set(zero_ecov_list) & set(snp_list_in_comb)))) - len(non_existant)       
    result=[score, covered, float(1-float(max(id_list))/score)]
    return result


def make_combinations(anchor_map, idd, comb_trshld, pval_param, z): 
    ### Get all DNAseI maps from table enhancer_maps EXCEPT FOR THE INPUT MAP 
    input_disease_id=idd
    zero_ecov_list=z
    cursor.execute("SELECT id, name FROM enhancer_maps e WHERE e.id!=%s AND NOT EXISTS (SELECT * FROM results WHERE disease_id=%s AND pval<%s and map_id=e.id) AND EXISTS (SELECT * FROM results WHERE map_id=e.id) ORDER BY id", (anchor_map, input_disease_id, pval_param))

    maps = cursor.fetchall()

    curr_map_id=0
    curr_map_id1=0
    curr_map_id2=0

    curr_map_name=""
    curr_map_name1=""
    curr_map_name2="" 

    query_maps=[]
    max_num=0
    max_maps=[] # map_ids that contain maximum number of requested elements (regions with ecov>0)
    all_ecov_maps=[]

    ### Identify all maps that could be potentially combined with first prioritized cadidate maps and also among each other put into all_ecov_maps 
    for m in maps:
        cursor.execute("SELECT h.snp_id, h.ecov FROM hits h WHERE h.disease_id=%s AND h.map_id=%s AND h.ecov>0 ORDER BY snp_id", (input_disease_id, m[0]))

        query_map=cursor.fetchall() 
        sum_nts=0   

        sums_dict={}    
        query_maps=[]
        for qm in query_map:
            query_maps.append(qm[0])
            sums_dict[qm[0]]=int(qm[1]) # snp coverage  (in a current disease)

        intersect=set(query_maps) & set(zero_ecov_list)
    
        if len(intersect) is not 0:
            all_ecov_maps.append(m[0])               

    tmp_list=[]
    combinations_dict={}
    i=0
    if comb_trshld > 2:
        for l in all_ecov_maps:
            for subset in itertools.combinations(all_ecov_maps, comb_trshld-1):
                if len(subset)>1:
                    tmp_list=list(subset)
                    tmp_list.append(anchor_map)
                    combinations_dict[i]=tmp_list
                    i+=1

    else:   ### combine all pairwise combinations with anchor map
        for l in all_ecov_maps:
            tmp_list=[]
            tmp_list.append(anchor_map)
            tmp_list.append(l)
            combinations_dict[i]=tmp_list
            i+=1    

    return combinations_dict


def merge_intervals(t_dict):
    merge_flag=0
    uniq_t_dict={}
    for k in t_dict.keys():
        if t_dict[k] and len(t_dict[k])>1:
            uniq_t_dict[k]=[]
            l=t_dict[k]

            # Sort tuples
            sorted_tuples = sorted(l)
            merged_list=[]
            tmp_el=0;
            cs, ce = sorted_tuples[0]
            ns, ne = sorted_tuples[1]
            j=0
            while j+1 < len(sorted_tuples):
                if ns <= ce:
                    merge_flag=1
                    merged_tuple = min(cs, ns), max(ce,ne)
                    cs,ce = merged_tuple    
                    j+=1
                    ns, ne = sorted_tuples[j]
                else:
                    merged_list.append((cs,ce))
                    cs, ce = (ns, ne)
            ns, ne = sorted_tuples[j]
            if ns <= ce:
                merged_tuple = min(cs,ns), max(ce,ne)
                merged_list.append(merged_tuple)
            else:
                merged_list.append((cs,ce))
                merged_list.append(sorted_tuples[j])

            if merged_list:
                uniq_t_dict[k] = merged_list
            else: 
                print "Error: Something went wrong with merging"
                sys.exit(1)
    return uniq_t_dict 


def sum_nts(t_dict):
    total_sum=0
    for k in t_dict.keys():
        tmp_list=t_dict[k]
        total_sum+= sum([i[1]-i[0] for i in tmp_list])
    
    return total_sum


def analyseCombinations(*args):
    if args[0].d is None:
        print "\nUse -h option to get more information."
        print "\nError: To few arguments."
        sys.exit(1) 
    else:
        input_disease_id=args[0].d

    # Set maximum number of maps to be combined
    if args[0].c is None:
        comb_trshld=2
    else:
        comb_trshld=args[0].c[0]
    
    # Set the number of processes that will run in parallel
    if args[0].t:
            ncpus=args[0].t[0]
            if ncpus > multiprocessing.cpu_count():
                print "WARNING: You don't seem to have enough processes available on your machine."
                print "Using "+str(multiprocessing.cpu_count()-1)+" processes instead of "+str(ncpus)+".\n"
    else:
            ncpus=multiprocessing.cpu_count()

    # Set the procentile of maps to be used as anchors
    if args[0].s:
        nt_param=args[0].s[0]
    else:
        nt_param=0.9

    # Number of recommeneded combinations of regulatory maps to print
    if args[0].n:
        n=args[0].n[0]
    else:
        n=3

    # Set a pval threshold so that only regulatory maps where input disease had pvalue > pval_param is used
    if args[0].p:
        pval_param=args[0].p[0]
    else:
        pval_param=0.001

    if args[0].m:
        improvement_coef=args[0].m[0]
    else:
        improvement_coef=5

    # --- End of input parsing --- #
    ################################

    ### Generate folder paths
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    DHS_maps_dir = os.path.join(parent_dir, "DHS_maps")
    tmp_dir = os.path.join(parent_dir, "tmp")
    combined_DHS_dir = os.path.join(tmp_dir, "combined_DHS_maps")

    # Create tmp/combined_DHS_dir folder if needed
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    if not os.path.exists(combined_DHS_dir):
        os.makedirs(combined_DHS_dir)

    i=0
    urr_map_name="0"
    curr_scores_dict={}
    stats_map_dict={}
    l=[]
    enhancer_maps_dict={}

   # remove any inconsistent data that might occured in previous runs..
    cursor.execute("SELECT id FROM enhancer_maps WHERE NOT EXISTS (SELECT count(*) FROM results WHERE map_id=id HAVING count(*)>1    )")
    to_remove=cursor.fetchall()
    for tr in to_remove:
        removeMap(tr[0])

    # Get names from all regulatory maps
    cursor.execute("SELECT id, name FROM enhancer_maps WHERE EXISTS (SELECT count(*) FROM results WHERE map_id=id HAVING count(*)>1) ORDER BY id")    # Previously tested regulatory maps stored in DB 
    emaps=cursor.fetchall()
    for e in emaps:
        enhancer_maps_dict[e[0]]=e[1]

    # Get total score and total number of nts from each map that did not achieve low pval
    cursor.execute("SELECT h.map_id, SUM(h.score), SUM(h.sizeb) FROM hits h WHERE h.disease_id=%s AND h.ecov>0 AND NOT EXISTS (SELECT * FROM results WHERE pval<%s AND disease_id=h.disease_id AND map_id=h.map_id) AND EXISTS (SELECT count(*) FROM results WHERE map_id=h.map_id HAVING count(*)>1) GROUP BY h.map_id", (input_disease_id, pval_param))

    stats_per_map=cursor.fetchall()
    for spm in stats_per_map:
        curr_scores_dict[int(spm[0])]=[int(spm[1]), int(spm[2])]
        l.append(int(spm[1]))
        #print int(spm[0]), curr_scores_dict[int(spm[0])]

    max_nt=max(l)
    nt_trshld=max_nt*nt_param

    ### Get the number of lSNP regions that have ecov=0 across ALL regulatory maps
    cursor.execute("SELECT snp_id FROM hits h WHERE h.disease_id=%s AND NOT EXISTS (SELECT * FROM hits WHERE snp_id=h.snp_id AND h.disease_id=disease_id AND ecov>0) AND map_id=1" % input_disease_id)
    non_existant_res=cursor.fetchall()
    non_existant=[]
    for ner in non_existant_res:
        non_existant.append(ner[0])

    max_results_list=[] 
    max_results_dict={}

    # Get the number of genomic regions for the input disease
    num_leaders=0
    cursor.execute("SELECT name, num_leaders FROM diseases WHERE id=%s" % input_disease_id)
    num_leaders=cursor.fetchall()[0]
    print "\nInput disease/trait: ", num_leaders[0], "("+str(input_disease_id)+")\n"

    zero_ecov_list=[]
    for spm in stats_per_map:
        if int(spm[1]) >= nt_trshld:
            stats_map_dict[spm[0]]=[int(spm[1]), int(spm[2])] 
            input_map_id=spm[0]

            # Get the map name
            cursor.execute("SELECT name FROM enhancer_maps WHERE id=%s" % input_map_id)
            input_map_name=cursor.fetchone()[0]

            # Get the number of genomic regions for the input disease
            num_leaders=0
            cursor.execute("SELECT name, num_leaders FROM diseases WHERE id=%s" % input_disease_id)
            num_leaders=cursor.fetchall()[0]
            print "Current seed map: ", input_map_name, "("+str(input_map_id)+")"

            # Get all regions from input map that have ecov=0 for the input (current) disease_id
            cursor.execute("SELECT snp_id, ecov, score FROM hits WHERE disease_id=%s AND map_id=%s AND ecov=0 ORDER BY snp_id", (input_disease_id, input_map_id))
            zero_ecov=cursor.fetchall()
            for ze in zero_ecov:    # in case additional info is needed (ecov and sizeb from the query), make zero_ecov_list to be a dict
                zero_ecov_list.append(ze[0])  

            adj_total_num=num_leaders[1] - len(non_existant)
            zero_ecov_list=list(set(zero_ecov_list) - set(non_existant))    

            if len(zero_ecov) is 0:
                print "Disease has all regions covered.. Exiting."

            combinations_dict={}
            prev_max_ncomb=0 

            combinations_dict=make_combinations(input_map_id, input_disease_id, comb_trshld, pval_param, zero_ecov_list)

            j=0
            scores_list=[]          # Keeping a list of #nts from all combinations in current anchor map
            covered_grs_dict={}     # Number of lSNP blocks from the current combination that are enriched with enhancers
            max_fold = []           # Combinations with max fold from each anchor
    
            # -------------------- #
            # Paralelize from here #
            # -------------------- #
            pool = multiprocessing.Pool(processes=ncpus)

            func = partial(analyseOneComb, combinations_dict, covered_grs_dict, curr_scores_dict, input_disease_id, num_leaders, zero_ecov_list, non_existant)
            try:
                results=pool.map(func, combinations_dict.keys())
            except:
                print "Pool.map error: Exiting..."
                sys.exit(1)
            pool.close()
            pool.join()
            k=0

            for i in results:
                scores_list.append(i[0])
                covered_grs_dict[k] = i[1]
                k+=1    
                max_fold.append(i[2])
                #print i[0], i[1], i[2]
            
            if max_fold:
                max_index=[m for m, r in enumerate(max_fold) if r == max(max_fold)]
                i=0
            
                # Get the max from all possible anchors
                max_results_list.append(round(max(max_fold),2))
                max_results_dict[spm[0]] = [combinations_dict[max_index[0]], round(max(max_fold),2)] #,scores_list[max_index[0]]]#

    # Remove combinations with same elements, if any
    del_list=[]
    for i in max_results_dict.keys():
        for j in max_results_dict.keys():
            if j>i:
                if len(set(max_results_dict[i][0]) & set(max_results_dict[j][0])) == comb_trshld:
                    del_list.append(j)

    for i in del_list:
        if i in max_results_dict.keys():
            del max_results_dict[i]

    scores_list1=[]

    for i in max_results_dict.keys():
        scores_list1.append(float(max_results_dict[i][1]))

    results_dict={}
    k=0
    min_old_p=1
    i=0

    for i in range(0,len(scores_list1)):
        max_index=[m for m, r in enumerate(scores_list1) if r == max(scores_list1)][0]
         
        ### Find combination with the hightest obtained score and ...
        for j in max_results_dict.keys():
            if max(scores_list1) == max_results_dict[j][1]:
                if k+1 <= n:
                    in_flag=1
                    f_tmp=""
                    name_str=""
                    min_old_p=1
                    c=0 

                    ### ... Concatenate two or more maps from selected combination of maps 
                    for el in max_results_dict[j][0]:
                        path_curr_map=os.path.join(DHS_maps_dir, enhancer_maps_dict[el]+"_input.csv")
                        if c==0:
                            # Create new file name by combining old ones
                            name_str=ntpath.basename(enhancer_maps_dict[el]).split('.')[0]
                            #print name_str
                            c=1
                        else:
                            name_str+='_'+ntpath.basename(enhancer_maps_dict[el]).split('.')[0]                         

                        # Find the smallest p value among combined regulatory maps --- this value is used later as an impact measurement
                        cursor.execute("SELECT pval FROM results WHERE map_id=%s AND disease_id=%s", (el, input_disease_id ))
                        get_pval=cursor.fetchone()[0]
                        min_old_p=get_pval if get_pval<min_old_p else min_old_p
                        with open(path_curr_map, "rb") as f1:
                            f_tmp+=f1.read()
                    with open(os.path.join(tmp_dir,"DEL.me"), "wb") as f:
                        f.write(f_tmp)

                    name_str+='_comb.csv'
                    print name_str, min_old_p # "combination id ", i, scores_list1[i]
                    ### Sort and merge overlapping intervals from combined maps
                    with open(os.path.join(tmp_dir,"DEL.me"), "rb") as f:
                        p1=Popen("sortBed -i", stdin=f, stdout=PIPE, shell=True)
                        p2=Popen("mergeBed -i", stdin=p1.stdout, stdout=PIPE, shell=True)
                        
                        with open(os.path.join(combined_DHS_dir, name_str), "wb") as fout:
                            fout.write(p2.stdout.read())
                    
                    test_path=os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_combination.py")

                    ### Start testing resulting map for a given disease/trait
                    with open(os.path.join(combined_DHS_dir, name_str), "rb") as f:
                        call('python '+test_path+' -i "%s" -d "%s"' % (str(os.path.join(combined_DHS_dir, name_str)),  str(input_disease_id)), shell=True)  

                    if os.path.exists(os.path.join(tmp_dir, "tmp.out")):
                        with open(os.path.join(tmp_dir, "tmp.out"), "a") as o:
                            o.write("\t"+str(min_old_p)+"\n")
                    else:
                        print "Temporary file does not exist, aborting."
                        sys.exit(1)
                    min_old_p=1
                    del max_results_dict[j]
                    k+=1

            scores_list1[max_index]=0

    ### Read results from temporary file
    s_line = []
    if os.path.exists(os.path.join(tmp_dir, "tmp.out")):
        with open(os.path.join(tmp_dir,"tmp.out"), "rb") as r:
            for line in r:
                s_line=line.split("\t")
                results_dict[s_line[0]]=[float(s_line[1]), float(s_line[2])]
    else:
        print "\nTesting was not performed, check your input parameters."
        sys.exit(1)
    print "\n"

    i=0 
    for result in results_dict.keys():
        min_old_p=results_dict[result][1]
        fold_change=float(min_old_p)/float(results_dict[result][0])
        if fold_change>=improvement_coef: #float(10):
            print str(results_dict[result][0])+" "+str(result)+" decreased p value %s times (old was %s)" % (int(fold_change), float(min_old_p))
        else:
            print str(results_dict[result][0])+" "+str(result)+" - p val did not improve (old was %s)" % float(min_old_p)
        i+=1

    ### Clean up a bit... 
    os.remove(os.path.join(tmp_dir, "tmp.out"))    
    os.remove(os.path.join(tmp_dir, "DEL.me"))
