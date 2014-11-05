#! /usr/bin/env python

########################################################################
# Extracts barcode and unique molecular index (UMI) from fastq record. #
# Stores in HDF5 file.
########################################################################
import sys
import os
import gzip
import bmh
import pdb
import signal
import sqlite3
import numpy as np
import Levenshtein as lev
from subprocess import Popen

#import bc_util
import pyximport; pyximport.install()
import bc_util_cy as bc_util
import fasta_cy as fasta

BARCODES = "/home/brad/lib/barcodes/bc2.txt"

cdef class extracter:
    cdef public char* fastq_fname
    cdef public char* dest_fname
    cdef int fastq_num_records
    #cdef int fastq_file
    cdef list barcodes
    cdef tuple bc_sub
    cdef tuple umi_sub

    def __init__(self, fastq_fname, dest_fname):

        ## FASTQ setup ##
        self.fastq_fname = fastq_fname
        print "Buffering fastq data..."

        # count fastq records if not already done so.
        # this is useful to predefine HDF5 table size
        if os.path.exists(fastq_fname + ".counts"):
            f = open(fastq_fname + ".counts")
            self.fastq_num_records = int(f.read())
        else:
            print "Counting fastq records..."
            self.fastq_num_records = fasta.count_records(fastq_fname)
            f = open(fastq_fname + ".counts", 'w')
            f.write(str(self.fastq_num_records))
            f.close()

        ## SQLite setup ##
        self.dest_fname = dest_fname
        #self.dest = 0

        # initialize dest if does not exist
        if not os.path.exists(self.dest_fname):
            initialize_dest(self.dest_fname)
        #self.conn = sqlite3.connect(self.dest_fname, isolation_level=None)
        #self.c = self.conn.cursor()

        ## Barcode setup ##
        barcodes_file = open(BARCODES)
        self.barcodes = []
        for l in barcodes_file:
            # expects format of name \t sequence
            self.barcodes.append(l.split()[1])

        ## Various other variables
        self.bc_sub = (5,11)
        self.umi_sub = (0,5)
        # self.i = 0
        # self.distance_min = 1E6
        # self.val = 0
        # self.ind = 0

    cpdef int read_fastq_sql(self):
        cdef int chunk_size = self.fastq_num_records / 4
        chunk_range = range(chunk_size)
        cdef int lineno = 0
        cdef int fastqno = 0
        cdef list name = ["" for x in chunk_range]
        cdef str seq = ""
        cdef str qual = ""
        cdef list bc = [0 for x in chunk_range]
        cdef list umi = ["" for x in chunk_range]

        fastq_file = fasta.open_gz(self.fastq_fname)
        conn = sqlite3.connect(self.dest_fname, isolation_level=None)
        c = conn.cursor()
        print "Reading " + self.fastq_fname

        cdef str l

        for l in fastq_file:
            #pdb.set_trace()
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
                self.insert_to_dest(name, bc, umi, c)
                fastqno = 0
            lineno += 1
        conn.commit()
        conn.close()

    # Search sequence [5:11] for barcodes
    # Returns index of barcode with minimum distance
    #cdef int get_barcode(self, seq, qual):
    cdef int get_barcode(self, seq, qual):
        #print seq, qual

        squal_min = fasta.min_qual(qual[self.bc_sub[0]:self.bc_sub[1]])
        #print squal_min
        # First, check qual. If minimum PHRED is less than 20 anywhere in barcode discard
        if squal_min <= 20: return -1

        sseq = seq[self.bc_sub[0]:self.bc_sub[1]]
        cdef int ind = 0
        #print sseq
        # Use Levenshtein distance to identify minimum scoring barcode
        ind = bc_util.min_barcode(self.barcodes, sseq, ind)
        return ind

    # Return candidate with best quality score
    # NOT NEEDED
    # def get_barcode_multiple(self, candidates, seq, qual):
    #     qual = [fasta.qual_to_int(q) for q in qual]
    #     qual_sum = [sum(q) for q in qual]
    #     return np.where(qual_sum == qual_sum.min())

    def get_umi(self, seq, qual):
        sub = (0, 5)
        min_qual = fasta.min_qual(qual[sub[0]:sub[1]])
        if (min_qual > 20):
            return seq[sub[0]:sub[1]]
        else:
            return 'N'

    cdef insert_to_dest(self, name, bc, umi, c):
        c.execute("BEGIN")
        for ind in xrange(len(name)):
            c.execute('''INSERT INTO data VALUES ('{0}', '{1}', '{2}', 0, 0)'''.format(name[ind], bc[ind], umi[ind]))
        c.execute("COMMIT")
        return 0


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
