import sys
#import sqlite3
import apsw

def printExplainQueryPlan(db, stmt, bindings):
    zExplain = "EXPLAIN QUERY PLAN " + stmt;

    ret = db.cursor().execute(zExplain, bindings)

    for line in ret.fetchall():
        print line
        #sys.stdout.write(line)

def read_sql(fname):
    with open(fname) as infile:
        return infile.read()
