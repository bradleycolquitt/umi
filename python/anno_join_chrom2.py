#! /usr/bin/env python

#############################################################################
#  Script to join barcode and UMI data from STRT data with gene annotations
#  Input database is created by load_to_db.py or load2db
#  Outputs tab-delimited text file with following columns:
#      barcode
#      aligned chromosome or reference
#      aligned read 5' position
#      total number of UMI at position
#      total number of unique UMI at position
#      gene ID
#      transcript ID
#      element
##############################################################################

import os
import sys
import sqlite3
import pdb
import argparse
import traceback as tb

from sql_utils import *
EXEC_PATH = os.path.dirname(os.path.realpath(__file__))

debug = True

def execute_join(db, build, statement, anno):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    anno = "/media/data/db/" + anno + ".db"
    c.execute("ATTACH DATABASE '{0}' AS anno".format(anno))

    c.execute("SELECT EXISTS(SELECT 1 FROM genes WHERE build=:build)", {"build":build})
    res = c.fetchone()[0]
    assert res == 1, "Given build is not present in annotation database."

    header = ""
    out = "_".join([db, build])
    create_table = ""
    sql = None
    if statement == "summary":
        sql = read_sql(EXEC_PATH + "/../sql/anno_join_summary.sql")

        header = "\t".join(["bc", "chrom", "lpos1", "lpos2", "rpos1", "rpos2",
                         "total_umi", "total_umi_same_strand", "unique_umi",
                         "gene_id", "transcript_id", "element"]) + "\n"
        out = out + "_summary.txt"
    elif statement == "full":
        #create_table = create_table_full

        sql = read_sql(EXEC_PATH + "/../sql/anno_join_full6.sql")
        #sql = sql % build

        header = "\t".join(["bc", "chrom", "lpos1", "lpos2", "rpos1", "rpos2", "umi",
                         "total_umi", "total_umi_same_strand",
                         "gene_id", "transcript_id", "element"]) + "\n"
        out = out + "_full.txt"

    dec = "y"
    # if os.path.exists(out):
    #     dec = raw_input("Output file exists. Overwrite [y/n]?")
    #     if not (dec == "y" or dec == "n"):
    #         print "Invalid input. [y/n]"
    #         sys.exit()



    #print build
    # Execute main join query
    if dec == "y":
        try:
            c.execute('DROP TABLE IF EXISTS %s' % build)
            #res = c.execute(create_table, {"build":build})
            create_sql = read_sql(EXEC_PATH + "/../sql/create_table_full2.sql") % build
            res = c.execute(create_sql)
            # Print out pla
            if (debug):
                printExplainQueryPlan(conn, sql, {"build":build})
            c.execute("PRAGMA synchronous = off")
            c.execute("PRAGMA journal_mode = MEMORY")
            res = c.execute(sql, {"build":build})
            # pdb.set_trace()
            # results = c.fetchall()

            # insertion_sql = "INSERT INTO {0} VALUES(?,?,?,?,?,?,?,?)".format(build)
            # c.execute("BEGIN TRANSACTION")
            # for r in results:

            #     c.execute(insertion_sql, r)
            # conn.commit()
        except sqlite3.OperationalError:
            print tb.print_exc()
            sys.exit(1)

        # Write to file
        # out = open(out, "w")
        # out.write(header)
        # for r1 in res.fetchall():
        #     out.write("\t".join(map(str, r1))+ "\n")
        # out.close()
        # conn.close()


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="db", help="Database containing umi, bc, and alignment data. As output by load_to_db.py")
    parser.add_argument(dest="build", help="Build from which to select annotation db elements.")
    parser.add_argument("-s", "--statement", dest="statement", choices=['full', 'summary'], default='summary' )
    parser.add_argument("-a", "--anno", default="anno", help="Annotation database prefix")
    args = parser.parse_args()

    execute_join(args.db, args.build, args.statement, args.anno)

if __name__ == "__main__":
    main(sys.argv)
