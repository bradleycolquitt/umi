#! /usr/bin/env python

#############################################################################
#  Script to join barcode and UMI data from STRT data with gene annotations
#  Input database is created by load_to_db.py
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

statement_summary = '''SELECT
                   data.bc,
                   reference.name,
                   align.position,
                   count(data.umi) as total_umi,
                   count(case when align.strand = anno.genes.strand then data.umi end)
                       as total_umi_same_strand,
                   count(distinct data.umi) as unique_umi,
                   anno.genes.gene_id,
                   anno.genes.transcript_id,
                   anno.genes.element
               FROM
                   data INNER JOIN align ON data.name = align.name
                        JOIN reference ON align.chrom = reference.chrom
                        JOIN anno.genes
                            ON reference.name = anno.genes.chrom
                                AND align.position BETWEEN anno.genes.start AND anno.genes.end
                                AND anno.genes.build=:build
               GROUP BY
                   data.bc, align.chrom, align.position

         '''

statement_full = '''SELECT
                   data.bc,
                   reference.name,
                   align.position,
                   data.umi,
                   count(data.umi) as umi,
                   count(case when align.strand = anno.genes.strand then data.umi end)
                       as total_umi_same_strand,
                   anno.genes.gene_id,
                   anno.genes.transcript_id,
                   anno.genes.element
               FROM
                   data INNER JOIN align ON data.name = align.name
                        JOIN reference ON align.chrom = reference.chrom
                        JOIN anno.genes
                            ON reference.name = anno.genes.chrom
                                AND align.position BETWEEN anno.genes.start AND anno.genes.end
                                AND anno.genes.build=:build
               GROUP BY
                   data.bc, data.umi, align.chrom, align.position

         '''

def execute_join(db, build, statement):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    c.execute("ATTACH DATABASE '/media/data/db/anno.db' AS anno")

    c.execute("SELECT EXISTS(SELECT 1 FROM genes WHERE build=:build)", {"build":build})
    res = c.fetchone()[0]
    assert res == 1, "Given build is not present in annotation database."

    header = ""
    out = ""
    if statement == "summary":
        statement = statement_summary
        header = "\t".join(["bc", "chrom", "position",
                         "total_umi", "total_umi_same_strand", "unique_umi",
                         "gene_id", "transcript_id", "element"]) + "\n"
        out = db + "_summary.txt"
    elif statement == "full":
        statement = statement_full
        header = "\t".join(["bc", "chrom", "position", "umi",
                         "total_umi", "total_umi_same_strand",
                         "gene_id", "transcript_id", "element"]) + "\n"
        out = db + "_full.txt"

    dec = "y"
    if os.path.exists(out):
        dec = raw_input("Output file exists. Overwrite [y/n]?")
        if not (dec == "y" or dec == "n"):
            print "Invalid input. [y/n]"
            sys.exit()
    # Execute main join query
    if dec == "y":
        try:
            res = c.execute(statement, {"build":build})
        except sqlite3.OperationalError:
            print tb.print_exc()

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
    args = parser.parse_args()

    execute_join(args.db, args.build, args.statement)

if __name__ == "__main__":
    main(sys.argv)
