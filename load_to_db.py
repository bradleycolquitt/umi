#! /usr/bin/env python

import sys
import os
import argparse

import pyximport; pyximport.install()
import bc_splitter_cy as bc
import bam_sql_cy as bs

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastq', dest='fastq', help='Original fastq file')
    parser.add_argument('-b', '--bam', dest='bam', help='Aligned reads. Must be indexed')

    args = parser.parse_args()

    db_path = "/media/data/strt_db/"
    db_fname = os.path.basename(args.fastq).split(".fastq")[0]
    db_fname = db_path + db_fname + ".db"

    try:
        ex = bc.extracter(args.fastq, db_fname)
        ex.read_fastq_sql()
    except:
        os.remove(db_fname)
        stop("Fastq load failed:" + tb.print_exc())

    try:
        bam_db = bs.bam_db(args.bam, db_fname)
        bam_db.fill_db()
    except:
        stop("BAM load failed:" + tb.print_exc())




if __name__ == "__main__":
    main(sys.argv)
