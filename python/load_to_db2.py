#! /usr/bin/env python

import sys
import os
import argparse
import traceback as tb

import pyximport; pyximport.install()
#import bc_splitter_cy as bc
import bam_sql_cy2 as bs

def main(argv):
    parser = argparse.ArgumentParser()
    #parser.add_argument('-f', '--fastq', dest='fastq', help='Original fastq file')
    parser.add_argument('-c', '--barcodes', dest='barcodes', help='Barcode file.')
    parser.add_argument('-b', '--bam', dest='bam', help='Aligned reads. Must be indexed')
    parser.add_argument('-o', '--output', help="Output database prefix")
    parser.add_argument("--nofilter_bc", default=False, action="store_true")
    parser.add_argument("--nofilter_umi", default=False, action="store_true")
    args = parser.parse_args()

    db_path = "/media/data/strt_db/"
    db_fname = ""
    if args.output:
        db_fname = args.output
    else:
        db_fname = os.path.basename(args.bam).split(".bam")[0]

    db_fname = db_path + db_fname + ".db"

    #try:
    #    ex = bc.extracter(args.fastq, db_fname, args.barcodes, args.nofilter_bc, args.nofilter_umi)
    #    ex.read_fastq_sql()
    #except:
    #    os.remove(db_fname)
    #    print "Fastq load failed:" + tb.print_exc()

    try:
        bam_db = bs.bam_db(args.bam, db_fname, args.barcodes, args.nofilter_bc, args.nofilter_umi)
        bam_db.fill_db()
    except:
        print "BAM load failed:" + tb.print_exc()




if __name__ == "__main__":
    main(sys.argv)
