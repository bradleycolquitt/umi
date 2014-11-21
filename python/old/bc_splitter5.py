#! /usr/bin/env python

########################################################################
# Extracts barcode and unique molecular index (UMI) from fastq record. #
# Stores in HDF5 file.
########################################################################
import sys
import os
import gzip
import bmh
#import fasta_cy as fasta
import fasta
import pdb
import signal
import sqlite3
import numpy as np
import Levenshtein as lev
from subprocess import Popen

import pyximport; pyximport.install()
from bc_util import *


BARCODES = "/home/brad/lib/barcodes/bc2.txt"

class extracter:
    def __init__(self, fastq_fname, dest_fname):
        try:
            ## FASTQ setup ##
            self.fastq_fname = fastq_fname
            print "Buffering fastq data..."

            self.fastq_file = fasta.open_gz(fastq_fname)

            # count fastq records if not already done so.
            # this is useful to predefine HDF5 table size
            # if os.path.exists(fastq_fname + ".counts"):
            #     f = open(fastq_fname + ".counts")
            #     self.fastq_num_records = int(f.read())
            # else:
            #     print "Counting fastq records..."
            #     self.fastq_num_records = fasta.count_records(fastq_fname)
            #     f = open(fastq_fname + ".counts", 'w')
            #     f.write(str(self.fastq_num_records))
            #     f.close()

            ## SQLite setup ##
            self.dest_fname = dest_fname
            self.dest = 0

            # initialize dest if does not exist
            if not os.path.exists(self.dest_fname):
                self.dest = initialize_dest(self.dest_fname)
            self.conn = sqlite3.connect(self.dest_fname, isolation_level=None)
            self.c = self.conn.cursor()

            ## Barcode setup ##
            self.barcodes_file = open(BARCODES)
            self.barcodes = []
            self.num_barcodes = 0
            for l in self.barcodes_file:
                # expects format of name \t sequence
                self.barcodes.append(l.split()[1])
                self.num_barcodes += 1

            ## Various other variables
            self.bc_sub = (5,11)
            self.umi_sub = (0,5)
            self.i = 0
            self.distance_min = 1E6
            self.val = 0
            self.ind = 0
            #self.bc_distances = np.zeros(self.num_barcodes)
        except AttributeError as e:
            print "AttributeError"
            #print "Attribute error({0}): {1}".format(e.errno, e.strerror)

    def read_fastq_sql(self):

        chunk_size = 100
        chunk_range = range(chunk_size)
        lineno = 0
        fastqno = 0
        name = ["" for x in chunk_range]
        seq = ""
        qual = ""
        bc = [0 for x in chunk_range]
        umi = ["" for x in chunk_range]

        print "Reading " + self.fastq_fname
        for l in self.fastq_file:
            if lineno % 4 == 0:
                name[fastqno] = l.split()[0].strip('@')
            if lineno % 4 == 1:
                seq = l
            if lineno % 4 == 3:
                qual = l
                bc[fastqno] = self.get_barcode(seq, qual)
                umi[fastqno] = self.get_umi(seq, qual)
                fastqno += 1
            if fastqno == chunk_size:
                self.insert_to_dest(name, bc, umi)
                fastqno = 0
            lineno += 1
        self.conn.commit()

    # Search sequence [5:11] for barcodes
    # Returns index of barcode with minimum distance
    def get_barcode(self, seq, qual):
        pdb.set_trace()
        squal_int = fasta.qual_to_int(qual[self.bc_sub[0]:self.bc_sub[1]])

        # First, check qual. If minimum PHRED is less than 20 anywhere in barcode discard
        if min(squal_int) <= 20: return 0

        sseq = seq[self.bc_sub[0]:self.bc_sub[1]]
        ind = 0

        # Use Levenshtein distance to identify minimum scoring barcode
        min_barcode(self.barcodes, sseq, ind)
        return ind

    # Return candidate with best quality score
    # NOT NEEDED
    # def get_barcode_multiple(self, candidates, seq, qual):
    #     qual = [fasta.qual_to_int(q) for q in qual]
    #     qual_sum = [sum(q) for q in qual]
    #     return np.where(qual_sum == qual_sum.min())

    def get_umi(self, seq, qual):
        sub = (0, 5)
        qual = fasta.qual_to_int(qual[sub[0]:sub[1]])
        if (min(qual) > 20):
            return seq[sub[0]:sub[1]]
        else:
            return 'N'

    def insert_to_dest(self, name, bc, umi):
        self.c.execute("BEGIN")
        for ind in xrange(len(name)):
            self.c.execute('''INSERT INTO data VALUES ('{0}', '{1}', '{2}', 0, 0)'''.format(name[ind], bc[ind], umi[ind]))
        self.c.execute("COMMIT")


# include functionality to read aligner info from MySQL db and add as attributes
#    * aligner
#    * alignment date
#    * command line
#    * barcode file
def initialize_dest(dest_fname):
    conn = sqlite3.connect(dest_fname)
    c = conn.cursor()

    c.execute('''CREATE TABLE data (name text, bc int, umi text, chrom int, position int)''')
    conn.close()
