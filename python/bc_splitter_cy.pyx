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

BARCODES = "/home/brad/lib/barcodes/bc2.txt"

cdef class extracter:
    cdef public char* fastq_fname
    cdef public char* dest_fname
    cdef int fastq_num_records
    cdef list barcodes
    cdef tuple bc_sub
    cdef tuple umi_sub
    cdef str sseq

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
        if not os.path.exists(self.dest_fname):
            self.initialize_dest()

        ## Barcode setup ##
        barcodes_file = open(BARCODES)
        self.barcodes = []
        for l in barcodes_file:
            self.barcodes.append(l.split()[1])  # expects format of name \t sequence

        ## Other variables ##
        self.bc_sub = (5,11)
        self.umi_sub = (0,5)
        self.sseq = "AAAAAA"
        #self.ind = 0

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
        squal_min = fasta.min_qual(qual[self.bc_sub[0]:self.bc_sub[1]])
        if squal_min <= 20: return -1

        self.sseq = seq[self.bc_sub[0]:self.bc_sub[1]]

        # Use Levenshtein distance to identify minimum scoring barcode
        return bc_util.min_barcode(self.barcodes, self.sseq)

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
