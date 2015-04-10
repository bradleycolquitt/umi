#! /usr/bin/env python

import os
import sys
import sqlite3
import pdb
import argparse
import traceback as tb

from sql_utils import *
EXEC_PATH = os.path.dirname(os.path.realpath(__file__))

debug = True

def execute_join(db, build):
    conn = sqlite3.connect(db)
    c = conn.cursor()

    header = ""
    out = "_".join([db, build])

    # header = "\t".join(["bc", "chrom", "lpos1", "lpos2", "rpos1", "rpos2",
    #                      "total_umi", "total_umi_same_strand",
    #                      "gene_id", "transcript_id", "gene_strand", "element"]) + "\n"
    header = "\t".join(["bc", "chrom", "lpos1", "lpos2", "rpos1", "rpos2", "strand",
                         "total_umi",
                         "gene_id", "transcript_id", "gene_strand", "element"]) + "\n"
    out = out + "_summary.txt"
    dec = "y"
    if os.path.exists(out):
        dec = raw_input("Output file exists. Overwrite [y/n]?")
        if not (dec == "y" or dec == "n"):
            print "Invalid input. [y/n]"
            sys.exit()

    # Execute main join query
    if dec == "y":
        try:
            sql = read_sql(EXEC_PATH + "/../sql/anno_join_full2.sql")
            if (debug):
                printExplainQueryPlan(conn, sql, {"build":build})
            res = c.execute(sql, {"build":build})

        except sqlite3.OperationalError:
            print tb.print_exc()
            sys.exit(1)

        #Write to file
        out = open(out, "w")
        out.write(header)
        for r1 in res.fetchall():
            out.write("\t".join(map(str, r1))+ "\n")
        out.close()
        conn.close()


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="db", help="Database containing umi, bc, and alignment data. As output by load_to_db.py")
    parser.add_argument(dest="build", help="Build from which to select annotation db elements.")

    args = parser.parse_args()

    execute_join(args.db, args.build)

if __name__ == "__main__":
    main(sys.argv)
