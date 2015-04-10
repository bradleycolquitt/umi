import sys
import sqlite3

def printExplainQueryPlan(db, stmt, bindings):
    zExplain = "EXPLAIN QUERY PLAN " + stmt;

    ret = db.execute(zExplain, bindings)

    for line in ret.fetchall():
        print line
        #sys.stdout.write(line)

def read_sql(fname):
    with open(fname) as infile:
        return infile.read()
