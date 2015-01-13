#! /usr/bin/env python

########################################################################
# Extracts barcode and unique molecular index (UMI) from fastq record. #
# Stores in SQLite database.
########################################################################
import sys
import os
import pdb
import sqlite3
import Levenshtein as lev

import pyximport; pyximport.install()
import bc_util_cy as bc_util
import fasta_cy as fasta

from cpython cimport bool

BARCODES = "/home/brad/lib/barcodes/"

# TODO
# filter out
cdef class extracter:
    cdef public char* fastq_fname
    cdef public char* dest_fname
    cdef int fastq_num_records
    cdef list barcodes
    cdef tuple bc_sub
    cdef tuple umi_sub
    cdef str sseq
    cdef bool nofilter_bc
    cdef bool nofilter_umi

    def __init__(self, fastq_fname, dest_fname, barcodes, nofilter_bc, nofilter_umi):

        ## FASTQ setup ##
        self.fastq_fname = fastq_fname

        ## SQLite setup ##
        self.dest_fname = dest_fname
        if not os.path.exists(self.dest_fname):
            self.initialize_dest()

        # count fastq records if not already done so.
        # this is useful to predefine sqlite chunk size
        if os.path.exists(fastq_fname + ".counts"):
            f = open(fastq_fname + ".counts")
            self.fastq_num_records = int(f.read())
        else:
            print "Counting fastq records..."
            self.fastq_num_records = fasta.count_records(fastq_fname)
            f = open(fastq_fname + ".counts", 'w')
            f.write(str(self.fastq_num_records))
            f.close()

        ## Barcode setup ##
        try:
            barcodes_file = open("".join([BARCODES, barcodes, ".txt"]))
        except IOError:
            raise("Barcode file not found.")
        self.barcodes = []
        for l in barcodes_file:
            self.barcodes.append(l.split()[1])  # expects format of name \t sequence

        ## Filter flags ##
        self.nofilter_bc = nofilter_bc
        self.nofilter_umi = nofilter_umi

        ## Other variables ##
        self.bc_sub = (5,11)
        self.umi_sub = (0,5)
        self.sseq = "AAAAAA"

    def initialize_dest(self):
        conn = sqlite3.connect(self.dest_fname)
        c = conn.cursor()
        c.execute('''CREATE TABLE data (name text, bc int, umi text, chrom int, position int)''')
        conn.close()

    ## Main function to read in fastq
    cpdef int read_fastq_sql(self):
        # Setup chunks for faster SQL insertion
        cdef int chunk_size = self.fastq_num_records / 4
        chunk_range = range(chunk_size)

        cdef int lineno = 0
        cdef int fastqno = 0
        cdef str seq = ""
        cdef str qual = ""

        cdef list name = ["" for x in chunk_range]
        cdef list bc = [0 for x in chunk_range]
        cdef list umi = ["" for x in chunk_range]

        print "Buffering fastq data..."
        fastq_file = fasta.open_gz(self.fastq_fname)
        conn = sqlite3.connect(self.dest_fname, isolation_level=None)
        c = conn.cursor()

        print "Reading " + self.fastq_fname
        cdef str l
        for l in fastq_file:
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

        print "Indexing..."
        c.execute('''
                  CREATE INDEX data_name ON data(name)
                  ''')
        c.execute('''
                  CREATE INDEX data_bc ON data(bc)
                  ''')
        c.execute('''
                  CREATE INDEX data_umi ON data(umi)
                  ''')
        conn.close()

    # Search sequence [5:11] for barcodes
    # Returns index of barcode with minimum distance
    cdef int get_barcode(self, seq, qual):
        # First, check qual. If minimum PHRED is less than 20 anywhere in barcode discard
        if not self.nofilter_bc:
            squal_min = fasta.min_qual(qual[self.bc_sub[0]:self.bc_sub[1]])
            if squal_min <= 20: return -1

        self.sseq = seq[self.bc_sub[0]:self.bc_sub[1]]

        # Use Levenshtein distance to identify minimum scoring barcode
        return bc_util.min_barcode(self.barcodes, self.sseq)

    def get_umi(self, seq, qual):
        sub = (0, 5)
        if not self.nofilter_umi:
            min_qual = fasta.min_qual(qual[sub[0]:sub[1]])
            if min_qual < 17:
                return 'N'
        return seq[sub[0]:sub[1]]

    cdef insert_to_dest(self, name, bc, umi, c):
        c.execute("BEGIN")
        for ind in xrange(len(name)):
            c.execute('''INSERT INTO data VALUES ('{0}', '{1}', '{2}', 0, 0)'''.format(name[ind], bc[ind], umi[ind]))
        c.execute("COMMIT")
