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

from sql_utils import printExplainQueryPlan

debug = True

statement_summary = '''
               SELECT
                   merge.bc,
                   anno.genes.chrom,
                   merge.hpos,
                   merge.tpos,
                   count(merge.umi) as total_umi,
                   count(case when merge.strand = anno.genes.strand then merge.umi end)
                       as total_umi_same_strand,
                   count(distinct merge.umi) as unique_umi,
                   anno.genes.gene_id,
                   anno.genes.transcript_id,
                   anno.genes.element
               FROM
                   merge, pos_rtree WHERE merge.rowid = pos_rtree.id
                   JOIN reference ON merge.tid = reference.tid
                   JOIN anno.genes ON reference.name = anno.genes.chrom
                        AND anno.genes.build=:build
                        AND
                        (
                        pos_rtree.lpos1 BETWEEN anno.genes.start AND anno.genes.end OR
                        pos_rtree.rpos2 BETWEEN anno.genes.start AND anno.genes.end
                        )
               GROUP BY
                   merge.bc, merge.tid, merge.rpos1

         '''
# FIX THIS!
statement_full = '''
               SELECT
                   merge.bc,
                   --anno.genes.chrom,
                   merge.lpos1,
                   merge.lpos2,
                   merge.rpos1,
                   merge.rpos2,
                   merge.umi,
                   count(merge.umi) as umi
                   --count(case when merge.strand = anno.genes.strand then merge.umi end)
                   --    as total_umi_same_strand,
                   --anno.genes.gene_id,
                   --anno.genes.transcript_id,
                   --anno.genes.element
               FROM
                  --merge
                  merge JOIN pos_rtree ON merge.rowid = pos_rtree.id
                       JOIN reference ON merge.tid = reference.tid
                       JOIN anno.genes ON reference.name = anno.genes.chrom
                           AND anno.genes.build=:build
                           AND
                           (
                           pos_rtree.lpos1 BETWEEN anno.genes.start AND anno.genes.end OR
                                pos_rtree.rpos2 BETWEEN anno.genes.start AND anno.genes.end
                           )

               GROUP BY
                   merge.bc, merge.umi, merge.tid, merge.lpos1
         '''

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
    if statement == "summary":
        statement = statement_summary
        header = "\t".join(["bc", "chrom", "lpos1", "lpos2", "rpos1", "rpos2",
                         "total_umi", "total_umi_same_strand", "unique_umi",
                         "gene_id", "transcript_id", "element"]) + "\n"
        out = out + "_summary.txt"
    elif statement == "full":
        statement = statement_full
        header = "\t".join(["bc", "chrom", "lpos1", "lpos2", "rpos1", "rpos2", "umi",
                         "total_umi", "total_umi_same_strand",
                         "gene_id", "transcript_id", "element"]) + "\n"
        out = out + "_full.txt"

    dec = "y"
    if os.path.exists(out):
        dec = raw_input("Output file exists. Overwrite [y/n]?")
        if not (dec == "y" or dec == "n"):
            print "Invalid input. [y/n]"
            sys.exit()

    # Print out plan
    if (debug):
        printExplainQueryPlan(conn, statement, {"build":build})

    # Execute main join query
    if dec == "y":
        try:
            res = c.execute(statement, {"build":build})
        except sqlite3.OperationalError:
            print tb.print_exc()
            sys.exit(1)

        # Write to file
        out = open(out, "w")
        out.write(header)
        for r1 in res.fetchall():
            out.write("\t".join(map(str, r1))+ "\n")

        conn.close()
        out.close()

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
