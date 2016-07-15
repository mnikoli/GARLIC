import MySQLdb as mdb
import sys

def mysql_connect():
    try:
        con=mdb.connect('localhost', 'root','milos', 'test_significance1')
        con.query('SET GLOBAL connect_timeout=172800')
        con.query('SET GLOBAL wait_timeout=172800')
        con.query('SET GLOBAL interactive_timeout=172800')
        return con
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit(1)

