#! /usr/bin/env python

#### TODO modify:
######   extract bc and umi from bam
######   use CIGAR info to filter and define

import os
import pysam
import pdb
import sqlite3

import pyximport; pyximport.install()
import bc_util_cy as bc_util

from cpython cimport array
from cpython cimport bool
#from array import array

BARCODES = "/home/brad/lib/barcodes/"

cdef class bam_db:
    cdef str bam_fname
    cdef object bam
    cdef tuple refs
    cdef long bam_counts
    cdef str dest_fname
    cdef object conn
    cdef object c
    cdef list barcodes
    cdef bool nofilter_bc
    cdef bool nofilter_umi

    def __init__(self, bam_fname, dest_fname, barcodes, nofilter_bc, nofilter_umi):
        self.bam_fname = bam_fname
        print self.bam_fname
        self.bam = pysam.AlignmentFile(bam_fname, 'rb')
        self.refs = self.bam.references
        self.bam_counts = self.bam.count()

        self.dest_fname = dest_fname

        self.conn = sqlite3.connect(self.dest_fname, isolation_level=None)
        self.c = self.conn.cursor()

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

    cpdef fill_db(self):
        # Insert BAM data into db
        print "Filling " + self.dest_fname
        self.c.execute('''CREATE TABLE IF NOT EXISTS align (name text,
                                                            chrom int,
                                                            position int,
                                                            strand int,
                                                            bc int,
                                                            umi text)''')
        cdef long readno = 0
        chunk_size = self.bam_counts / 10
        cdef object read
        cdef list cigar
        cdef int pos
        cdef str seq
        cdef array.array qual
        cdef char[:] cqual
        cdef int bc
        cdef str umi

        self.c.execute("BEGIN TRANSACTION")

        for read in self.bam.fetch():
            if not read.is_unmapped:
                # CIGAR>25 indicator of poorly structured barcode sequence
                cigar = read.cigartuples
                if cigar[0][1] > 25: continue

                if read.is_reverse:
                    pos = read.alen
                else:
                    pos = read.pos
                seq = read.query_sequence
                qual = read.query_qualities
                cqual = qual
                bc = bc_util.get_barcode(seq,
                                         cqual,
                                         (6,11),
                                         self.nofilter_bc,
                                         self.barcodes)
                umi = bc_util.get_umi(seq,
                                      cqual,
                                      (0,5),
                                      self.nofilter_umi)

                self.c.execute('''INSERT INTO align VALUES (?, ?, ?, ?, ?, ?)''', (read.qname,
                                                                             read.rname,
                                                                             pos,
                                                                             read.is_reverse,
                                                                             bc,
                                                                             umi))
                readno += 1
            if readno == chunk_size:
                self.c.execute("COMMIT")
                self.c.execute("BEGIN TRANSACTION")
                readno = 0
        self.c.execute("COMMIT")

        # Create index for name, chrom, position
        print "Indexing..."
        self.create_indices()

        self.create_reftable()

    def create_reftable(self):
        # Create references table containing key for refid and human-readable chromosome label
        print "Create Ref id table"
        self.c.execute("BEGIN TRANSACTION")
        self.c.execute('''CREATE TABLE IF NOT EXISTS reference (name text, chrom int)''')
        for idx, ref in enumerate(self.refs):
            self.c.execute('''INSERT INTO reference VALUES ('{0}', {1})'''.format(ref, idx))
        self.c.execute("COMMIT")

    def create_indices(self):
        self.c.execute('''
                       CREATE INDEX align_name ON align(name)
                       ''')
        self.c.execute('''
                       CREATE INDEX align_chrom_position ON align(chrom, position)
                       ''')
