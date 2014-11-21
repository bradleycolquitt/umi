#! /usr/bin/env python

## General pipeline for processing BC/UMI sequences

import sys
import os
import pdb
import subprocess
import sqlite3

import pyximport; pyximport.install()
import bc_splitter_cy as bc
import bam_sql_cy as bs


def main(argv):
    ROOT="/home/brad/src/seq_utils/python/test/"
    BASENAME = argv[1]
    FASTQ=ROOT + BASENAME + ".fastq.gz"
    DB=ROOT + "../../db/" + BASENAME + "_bc.db"
    BAM=FASTQ.split(".")[0]+".bam"
    #pdb.set_trace()
    ## Align reads to genome using bowtie2
    if not os.path.exists(BAM):
        print "Aligning..."
        cmd_args = "bowtie2 -p 10 -x taeGut1 -U {0} -S {1}; process_sam.sh {1}".format(FASTQ, FASTQ.split(".")[0]+".sam")
        print cmd_args
        p = subprocess.Popen(cmd_args, shell=True)
        p.wait()

    ## Load BC and UMI from fastq into HDF5
    ex = bc.extracter(FASTQ, DB)
    ret = ex.read_fastq_sql()
    #ex.conn.close()

    ## Update HDF5 with read aligned chrom and pos using read name as key
    bs1 = bs.bam_db(BAM, DB)
    print "fill_db"
    bs1.fill_dest()
    #bs1.conn.close()

    print "join"
    ba = bs.bc_align(ROOT + "../../db/" + BASENAME)
    ba.join_bc_align()

if __name__ == "__main__":
    main(sys.argv)
