#!/usr/bin/python

import mySQLConnect
from mySQLConnect import mysql_connect
# Connect to GARLIC DB
con=mysql_connect()
cursor=con.cursor()

cursor.execute("SELECT COUNT(*) FROM diseases WHERE name=%s" % disease_name);

