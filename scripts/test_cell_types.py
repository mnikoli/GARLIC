#!/usr/bin/python 
########################################################
# Script that prints tested regulatory maps ordered by 
# predicted etiological connection with a given disease                                   
########################################################

import sys
import os
import xlwt
import MySQLdb
from mySQLConnect import mysql_connect
#import MySQLdb as mdb
    
con=mysql_connect()
cursor=con.cursor()

def testCellTypes(*args):

    if args[0].did is None:
        print "Input error: Missing disease id"
        print "Usage: garlic testCellTypes <disease_id>"
        sys.exit(1)
    else:
        curr_disease_id=args[0].did
    
    print curr_disease_id, int(curr_disease_id), args[0].did

    ### Get disease/trait details
    try:
        cursor.execute("SELECT name, id, num_leaders, num_followers FROM diseases WHERE id=%s" % curr_disease_id)
    except:
        print "Error: Provided disease_id does not exist in database."
        print "\n Try './garlic viewTableData -n diseases' to get the list of available disease ids"
        sys.exit(1)

    info=cursor.fetchall()[0]

    info_row = "# "+str(info[0])+", id "+str(info[1])+", lSNPs "+str(info[2])+", fSNP "+str(info[3])
    print info_row
    init_num_set=[]
    generated_leader_ids=[]
    
    results_xls = xlwt.Workbook()
    sheet1 = "list"
    all_results = results_xls.add_sheet(sheet1)
    	
    rows=0
    cols=0

    ### NOTE This part is hard-coded
    if info[2]<7:
        print  "Results not available for '"+str(info[0])+"', as it has less than 7 GRRs."
        print "\n Try './garlic viewTableData -n diseases' to get the list of available disease ids"
        sys.exit(1)

    ### Get all trait details from table results
    cursor.execute("SELECT name, pval, disease_score, count, corrected_cell_pval FROM results, enhancer_maps WHERE disease_id=%s AND map_id=id ORDER BY pval" % curr_disease_id)	
    rows=cursor.fetchall()
    count=len(rows)
    
    ### Create new Excel sheet
    results_xls = xlwt.Workbook()
    sheet1 = str("list")
    all_results = results_xls.add_sheet(sheet1)
    
    col_names = ["Cell_Type", "Corrected_P_Value", "P_Value", "Disease_Score", "Count"]
    
    # Print trait details to XLS
    all_results.write(0,0,info_row)
    
    # Print column names
    for cols in range(0, len(col_names)):
        all_results.write(2, cols, col_names[cols])
    xrows=3
    p=0
    xcols=0
    
    col_results=[]
    for row in rows:
        xcols=0
        col_results=[str(row[0]), row[4], str(row[1]), str(row[2]),str(row[3])]
        for cols in range(0,len(col_names)):
            all_results.write(xrows, cols, col_results[xcols])
            xcols+=1
    
        xrows+=1
        p+=1
    
    parent_dir=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    map_dir = os.path.join(os.path.join(parent_dir, "reports"), "testedCellTypes")
    report_path = os.path.join(map_dir, str(info[0])+".xls")
    
    results_xls.save(report_path)
    print "Report has been generated in: "+report_path
