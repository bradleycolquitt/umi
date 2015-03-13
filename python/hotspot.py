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

# pos_select = '''
# SELECT
#     align.bc,
#     reference.name,
#     align.position,
#     count(align.umi) as total_umi
# FROM
#     align JOIN ON align.tid = reference.tid
# GROUP BY
#     align.tid, align.position, align.bc
# WHERE
#     reference.name = positions.ref_name AND
#     align.position = positions.pos
# (
# SELECT *
# FROM
#     (
#     SELECT
#         reference.name,
#         align.position,
#         count(align.umi) as total_umi
#     FROM
#         align JOIN reference ON align.tid = reference.tid
#     GROUP BY
#         align.tid, align.position
#     ) AS inner_table
#     WHERE
#         total_umi > 10
# ) AS positions
# '''

pos_select = '''
SELECT *
FROM
    (
    SELECT
        align.tid,
        align.position,
        count(align.umi) as total_umi
    FROM
        align JOIN reference ON align.tid = reference.tid
    GROUP BY
        align.tid, align.position
    ) AS inner_table
    WHERE
        total_umi > (?)
'''
pos_select_summary = '''
    SELECT
        align.bc,
        reference.name,
        align.position,
        count(*)
    FROM
        align
    WHERE
        align.tid = (?) AND align.position = (?)
    GROUP BY
        align.tid, align.position, align.bc

'''
pos_select_full = '''
    SELECT
        align.bc,
        reference.name,
        align.position,
        align.umi,
        count(*)
    FROM
        align JOIN reference ON align.tid = reference.tid
    WHERE
        align.tid = (?) AND align.position = (?)
    GROUP BY
        align.tid, align.position, align.bc, align.umi
'''

def execute_join(db, thresh, statement):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    header = ""
    out = db + "_hotspots"
    if statement == "summary":
        pos_select2 = pos_select_summary
        header = "\t".join(["bc", "chrom", "position", "count", "group"]) + "\n"
        out = out + "_summary.txt"
    elif statement == "full":
        pos_select2 = pos_select_full
        header = "\t".join(["bc", "chrom", "position", "umi", "umi_count", "group"]) + "\n"
        out = out + "_full.txt"

    dec = "y"
    if os.path.exists(out):
        dec = raw_input("Output file exists. Overwrite [y/n]?")
        if not (dec == "y" or dec == "n"):
            print "Invalid input. [y/n]"
            sys.exit()

    # Execute main join query
    if dec == "y":

        out = open(out, "w")
        out.write(header)

        max_dist = 100
        curr_value = 0
        group = 0
        try:
            results = c.execute(pos_select, (thresh,))
            for result in results.fetchall():
                results2 = c.execute(pos_select2, (result[0], result[1]))
                for res in results2.fetchall():

                    if (res[2] - curr_value) > max_dist:
                        group += 1
                    curr_value = res[2]

                    out.write("\t".join(map(str, res) + [str(group)])+ "\n")
        except sqlite3.OperationalError:
            print tb.print_exc()
            sys.exit()
        except:
            print tb.print_exc()

        conn.close()
        out.close()

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="db", help="Database containing umi, bc, and alignment data. As output by load_to_db.py")
    parser.add_argument("-t", dest="thresh", type=int)
    parser.add_argument("-s", "--statement", dest="statement", choices=['full', 'summary'], default='summary' )
    args = parser.parse_args()

    execute_join(args.db, args.thresh, args.statement)

if __name__ == "__main__":
    main(sys.argv)
