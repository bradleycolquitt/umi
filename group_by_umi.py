#! /usr/bin/env python

# group by chrom
# group by position
# count total number of umi
# count number of unique umi
# join position and annotations.genes

import sys
import sqlite3
import pdb

st1 = '''SELECT
             data.bc, align.chrom, align.position, count(data.umi) as total_umi, count(distinct data.umi) as unique_umi
         FROM
             data INNER JOIN align WHERE data.name = align.name JOIN
         GROUP BY
             data.bc, align.chrom, align.position
         '''

st2 = '''SELECT
             data.bc, reference.name, align.position, count(data.umi) as total_umi, count(distinct data.umi) as unique_umi, anno.genes.transcript_id
         FROM
             data INNER JOIN align ON data.name = align.name
                  JOIN reference ON align.chrom = reference.chrom
                  INNER JOIN anno.genes
                      ON reference.name = anno.genes.chrom
                          AND align.position BETWEEN anno.genes.start AND anno.genes.end
         GROUP BY
             data.bc, align.chrom, align.position
         '''

def execute_join(db):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    c.execute("ATTACH DATABASE '/media/data/db/anno.db' AS anno")
    res = c.execute(st2)

    out = open(db + "_test.txt", "w")
    out.write("\t".join(["bc", "chrom", "position", "total_umi", "unique_umi", "transcript_id"]) + "\n")
    for r1 in res.fetchall():
        out.write("\t".join(map(str, r1))+ "\n")
        #pass

    conn.close()
    out.close()

def main(argv):
    execute_join(argv[1])
    #pdb.set_trace()


if __name__ == "__main__":
    main(sys.argv)
