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

import sys
import sqlite3
import pdb
import traceback as tb

statement = '''SELECT
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

def execute_join(db, build):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    c.execute("ATTACH DATABASE '/media/data/db/anno.db' AS anno")

    c.execute("SELECT EXISTS(SELECT 1 FROM genes WHERE build=:build)", {"build":build})
    res = c.fetchone()[0]
    assert res == 1, "Given build is not present in annotation database."

    # Execute main join query
    try:
        res = c.execute(statement, {"build":build})
    except sqlite3.OperationalError:
        print tb.print_exc()

    # Write to file
    out = open(db + "_test.txt", "w")
    out.write("\t".join(["bc", "chrom", "position",
                         "total_umi", "total_umi_same_strand", "unique_umi",
                         "gene_id", "transcript_id", "element"]) + "\n")
    for r1 in res.fetchall():
        out.write("\t".join(map(str, r1))+ "\n")

    conn.close()
    out.close()

def main(argv):
    execute_join(argv[1], argv[2])

if __name__ == "__main__":
    main(sys.argv)
