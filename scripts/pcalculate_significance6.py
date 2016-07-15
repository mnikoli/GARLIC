#!/usr/bin/python 
###################################################################
# Parallelized script that performs the random sampling procedure #
###################################################################
import sys
import csv
import random
import threading
import multiprocessing
from multiprocessing import Process, Value
import os
import xlwt
from mySQLConnect import mysql_connect
import MySQLdb as mdb

from rpy2.robjects.packages import importr
stats = importr('stats', robject_translations={'format_perc': '_format_perc'})    
from rpy2.robjects.vectors import FloatVector

from decimal import Decimal
from scipy import stats as s
import numpy

def calculate_disease_score(leader_subset, leader_dict, curr_disease_id):
    sum_sizeb=0
    sum_ecov=0
    ecov=0
    disease_score=0
    scores=[]
    scores_dict={}
    num_nt=0
    num_nt_ecov=0
    
    ## Get the enhancer coverage score from all LD regions
    for ls in leader_subset:    
        scores=[int(leader_dict[ls][0]), int(leader_dict[ls][5])]
        scores_dict[ls]=scores

    score_sort=sorted(scores_dict.items(), key=lambda(k,v):v[0])

    ### Trimmed mean is used as a score (20% is trimmed from each end of the list, then mean is calculated)
    scores_keys=[]
    for ss in range(0, len(score_sort)):
        scores_keys.append(score_sort[ss][0])

    trim_par=0.2
    trim_size=int(trim_par*len(scores_keys))

    del scores_keys[len(scores_keys)-trim_size:]
    del scores_keys[:trim_size]

    trimmed_scores_keys=scores_keys

    ### Calculate the average number of nts overlapping regulatory regions  and the total number of overlapping nts after trimming 
    num_nt=[]
    num_nt_ecov=[]
    for k in trimmed_scores_keys:
        num_nt_ecov.append(scores_dict[k][0])
        num_nt.append(scores_dict[k][1])

    avg=numpy.average(num_nt_ecov)
    disease_score=float(avg)/sum(num_nt)

    return disease_score

def worker(d, sem, leader_dict, curr_map_id, init_repeats, num_bins):
    p = multiprocessing.Process(target=calculate_p, args=(d, sem, leader_dict, curr_map_id, init_repeats, num_bins))
    p.start()
    p.join()
    sem.release()

def calculate_p(d, sem, leader_dict, curr_map_id, init_repeats, num_bins):
    repeats=init_repeats
    curr_disease_id=d[0]
    curr_disease_name=d[1]
    leaders_in_disease=d[2]
    curr_num_followers=d[3]

    ### Set seed for the current disease by using curr_disease_id
    random.seed(curr_disease_id)
    con=mysql_connect()
    cursor=con.cursor()

    ### Initial num set gets first ids from all snps..
    init_num_set=[]
    init_num_set=leader_dict.keys()

    x=len(init_num_set)
    disease_num_set=[]
    rest=[]

    ### Get ids from all leaders that are associated with the current disease (from > v3.0 - they will be excluded in the random sampling procedure)
    try:
        cursor.execute("SELECT gr_id, snp_id FROM genomic_regions where disease_id=%s" % curr_disease_id)
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)

    excluded_grs=cursor.fetchall()
    exclude={}
    
    for e in excluded_grs:
        exclude[int(e[1])]=int(e[0])    #snps and grs
        disease_num_set.append(int(e[0])) #grs
        init_num_set.remove(int(e[0]))
#print len(leader_dict), int(e[0])

    # Calculate observed score for the current disease
    disease_score=calculate_disease_score(disease_num_set, leader_dict, curr_disease_id)
    print curr_disease_id, curr_disease_name, disease_score

    # Eliminate duplicate regions (same startb and endb) represented by the same lSNP but across different diseases
    # This set is different for each disease, that why it has to be recalculated every time
    dup_gr_snps_dict={}
    cursor.execute("SELECT DISTINCT snp_id FROM genomic_regions gr WHERE EXISTS (SELECT * FROM genomic_regions WHERE gr.snp_id=snp_id AND gr.disease_id!=disease_id AND gr.startb=startb AND gr.endb=endb) order by snp_id")
    dup = cursor.fetchall()
    #print "Num. of lSNPs with exact duplicate regions ", len(dup)

    # lSNP ids are used as dict keys
    tmp_dict={}
    for d in dup:
        tmp_dict[int(d[0])]=[]

    cursor.execute("SELECT snp_id, gr_id FROM genomic_regions gr WHERE EXISTS (SELECT * FROM genomic_regions WHERE gr.snp_id=snp_id AND gr.disease_id!=disease_id AND gr.startb=startb AND gr.endb=endb) order by snp_id")

    dup = cursor.fetchall()
    for d in dup:
        tmp_list=tmp_dict[int(d[0])]    
        tmp_list.append(int(d[1]))
        tmp_dict[int(d[0])]=tmp_list
    
    # tmp_dict[lsnp]=[gr1, gr2,..,grn]
    # go through snp ids...
    for k in tmp_dict.keys():
        tmp_list=tmp_dict[k]    # grab the list of grs
        if k in exclude.keys() and exclude[k] in tmp_list:    # if it's a snp from the current disease
            for t in tmp_list:
                if t in init_num_set:
                    init_num_set.remove(t)
        else:
            tmp_list.pop()
            for t in tmp_list:
                if t in init_num_set:
                    init_num_set.remove(t)

    leaders_in_bins={}    # dict with snp_ids in separate lists (bin_id is used as a dict key) 

    ### Initialize dict elements (bins) as lists
    for i in range(1,num_bins+1):    
        leaders_in_bins[i]=[]
    
    # Distribute SNPs across the bins according to their GR size
    for i in init_num_set:
        #print i,", ", unique_leader_grs_dict[i][1]
        leaders_in_bins[leader_dict[i][1]].append(i)            

    bins_list=[]
    bins_cc=[]    # content count    from bins
    for dns in disease_num_set:
        bins_list.append(leader_dict[dns][1])
    for i in range(1,num_bins+1):
        bins_cc.append(bins_list.count(i)) # this list contains number of elements in each of the #bin_size bins    
    
    if disease_score != 0:
        j=1
        k=0
        e_flag=0

        ### Randomly draw a certain number of leaders for many many times in order to estimate significance
        for j in range(1, repeats+1):
            sampled=[]    
            sim_score=0
            for i in range(1,num_bins+1):
                sampled_tmp=[]
                sampled_tmp=random.sample(leaders_in_bins[i], bins_cc[i-1])
                sampled=sampled+sampled_tmp

            sim_score=calculate_disease_score(sampled, leader_dict, curr_disease_id)
            if sim_score >= disease_score:
                k+=1
        try:
            cursor.execute("INSERT INTO results (disease_id, pval, disease_score, count, map_id) VALUES(%s,%s,%s,%s,%s)", (curr_disease_id, round(float(k+1)/(repeats+1), 10), round(disease_score,10), k, curr_map_id))
            con.commit()
        except mdb.Error, e:
            if con:
                 con.rollback()
                 print "Error %d: %s" % (e.args[0], e.args[1]) 
                 sys.exit(1)
    else:
        try:
            cursor.execute("INSERT INTO results (disease_id, pval, disease_score, count, map_id) VALUES(%s,%s,%s,%s,%s)", (curr_disease_id, 1, round(disease_score,10), -1, curr_map_id))
            con.commit()
        except mdb.Error, e:
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.exit(1)

def main_pcalculate_significance(map_id, name, ltrshld, ncpus, init_repeats):
    curr_map_id=map_id
    con=mysql_connect()
    cursor=con.cursor()
    
    init_num_set=[]
    generated_leader_ids=[]

    ###########################
    # INITIAL PARAMETERS ! ! !#
    ###########################

    num_bins=5 # do not change this
    
    print "Leaders per disease treshold: ", ltrshld
    print "Number of threads: ", ncpus
    print "Number of iterations: ", init_repeats, "\n"        

    ### Get the parent dir
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    ### Check for results dir 
    results_dir = os.path.join(os.path.join(parent_dir, "reports"), "testedDiseases")
    if not os.path.exists(results_dir):
        print "Error: Directory "+results_dir+" does not exist!"
        sys.exit(1)

    ### Check for results/id_name dir
#    map_dir = os.path.join(results_dir, str(curr_map_id)+"_"+name)
#    if not os.path.exists(map_dir):
#        print "Error: Directory "+map_dir+" does not exist!"
#        sys.exit(1)

    curr_follower_name="0"
    curr_leader_name="0"

    all_leaders_dict={}
    leader_dict={}

    ### Get ids from all leaders with additional info (ecov, bin_id, chr_id, startb, endb - I added last 3 terms)        # 1-bin, 2-chr_id, 3-startb, 4-endb, 5-sizeb, 6-snp_id
    try:
        cursor.execute("SELECT s.id, h.ecov, gr.bin, l.chr_id, gr.startb, gr.endb, gr.sizeb, gr.gr_id  FROM snps s, hits h, leaders l, genomic_regions gr WHERE s.leader=%s AND s.id=h.snp_id AND map_id=%s AND l.snp_id=s.id AND gr.snp_id=s.id AND h.disease_id=gr.disease_id", (1, curr_map_id));
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)
    fetched_leaders=cursor.fetchall()

    for fl in fetched_leaders:
        tmp_list=[]
        tmp_list=[fl[1], int(fl[2]), int(fl[3]), str(fl[4]), str(fl[5]), str(fl[6]), int(fl[0])]
        all_leaders_dict[int(fl[7])]=tmp_list

    leader_dict={}                        # NOTE: tmp fix
    leader_dict=all_leaders_dict
    
    # Get the list of all diseases that pass the given cutoff (#lSNPs)
    try:
        cursor.execute("SELECT id, name, num_leaders, num_followers FROM diseases WHERE num_leaders>=%s ORDER BY id" % ltrshld)
    except mdb.Error, e:
                print "Error %d: %s" % (e.args[0], e.args[1])
                sys.exit(1)
    diseases=cursor.fetchall()

    #-----------------------#
    # Parallelize from here #
    #-----------------------#

    k=Value('i',0)

    workers = list()
    sem = multiprocessing.Semaphore(ncpus)
    for d in diseases:
        workers.append(threading.Thread(target=worker, args=(d, sem, leader_dict, curr_map_id, init_repeats, num_bins)))
        sem.acquire()
        workers[-1].start()

    for w in workers: w.join()

#    if not status_file:
#        print "Error: Cannot open file "+status_path
#        sys.exit(1)
#    else:
#        status_file.write(str(main_pcalculate_significance.__name__)+": Successfully calculated disease significance for a map "+str(name)+")\n")

    ### Commit from main so that results are visible 
    con.commit()

    results_xls = xlwt.Workbook()
    sheet1 = "All results"
    all_results = results_xls.add_sheet(sheet1)
    
    rows=0
    cols=0

    col_names = ["Disease_ID", "Disease_Name", "Corrected_P_Value", "P_Value", "Disease_Score", "Num_Leaders", "Num_Followers", "Count"]
    
    ### Print column names in sheet1
    for cols in range(0,len(col_names)):
        all_results.write(0, cols, col_names[cols])

    cursor.execute("SELECT COUNT(*) FROM results WHERE map_id=%s" % curr_map_id)
    count=cursor.fetchone()[0]    

    cursor.execute("SELECT r.disease_id, d.name, r.pval, r.disease_score, d.num_leaders, d.num_followers, r.count FROM results r, diseases d WHERE r.disease_id=d.id AND map_id=%s ORDER BY pval" % curr_map_id)    

    results=cursor.fetchall()
    col_results=[]
    p_list=[]
    
    for result in results:
        p_list.append(result[2])

    p_adjust=stats.p_adjust(FloatVector(p_list), method='BH')
    
    rows=1
    p=0
    for result in results:
        cols=0
        results_tmp1 = [str(result[2]), str(1)][result[2] == 1.0]
        results_tmp0 = [str(result[3]), str(0)][result[3] == 0.0]
        col_results=[str(result[0]), str(result[1]),p_adjust[p], results_tmp1, results_tmp0, str(result[4]), str(result[5]), str(result[6])]
        cursor.execute("UPDATE results SET corrected_pval=%s WHERE map_id=%s AND disease_id=%s", (round(p_adjust[p],10), curr_map_id, str(result[0])))
        for cols in range(0,len(col_names)):
            all_results.write(rows, cols, col_results[cols])
            cols+=1    
        rows+=1
        p+=1
    con.commit()
    file_path=os.path.join(results_dir, name+"_"+str(curr_map_id)+".xls")
    results_xls.save(file_path)    
    print "Report has been generated in: "+file_path
