#! /usr/bin/env python

import os
import sys
import pysam
import pdb
import sqlite3
import re

import pyximport; pyximport.install()
import bc_util_cy as bc_util

from cpython cimport array
from cpython cimport bool

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
                                                            instrument text,
                                                            flowcell text,
                                                            chrom int,
                                                            position int,
                                                            strand int,
                                                            bc int,
                                                            umi text)''')
        cdef long readno = 0
        chunk_size = self.bam_counts / 10
        cdef object read
        cdef list cigar
        cdef str cigarstring
        cdef list tags
        cdef str qname
        cdef list qname_list
        cdef str instrument
        cdef str flowcell
        cdef int pos
        cdef str seq
        cdef array.array qual
        cdef char[:] cqual
        cdef int bc
        cdef str umi

        self.c.execute("BEGIN TRANSACTION")

        re_cigar = re.compile("[IDN]")

        for read in self.bam.fetch():
            if not read.is_unmapped:
#                if readno == 10: sys.exit()
                # Filter out reads with insertions, deletions, or reference skips
                cigarstring = read.cigarstring
                if re_cigar.search(cigarstring): continue

                # Filter out multimappers
                tags = read.tags
                if tags[0][1] > 1: continue

                qname = read.qname
                qname_list = qname.split(":")
                instrument = qname_list[0]
                flowcell = qname_list[2]

                if read.is_reverse:
                    # reference_end points to one past last aligned residue.
                    # Will convert to last residue with switch to 1-based indexing
                    pos = read.reference_end
                else:
                    pos = read.reference_start + 1

                seq = read.query_sequence
                qual = read.query_qualities
                cqual = qual

                bc = bc_util.get_barcode(seq,
                                         cqual,
                                         (5,11),
                                         self.nofilter_bc,
                                         self.barcodes)
#                print bc
                umi = bc_util.get_umi(seq,
                                      cqual,
                                      (0,5),
                                      self.nofilter_umi)

                self.c.execute('''INSERT INTO align VALUES (?, ?, ?, ?, ?, ?, ?, ?)''', (read.qname,
                                                                                         instrument,
                                                                                         flowcell,
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
        self.create_indices()

        # Create references table containing key for refid and human-readable chromosome label
        self.create_reftable()

    def create_reftable(self):

        print "Create Ref id table"
        self.c.execute("BEGIN TRANSACTION")
        self.c.execute('''CREATE TABLE IF NOT EXISTS reference (name text, chrom int)''')
        for idx, ref in enumerate(self.refs):
            self.c.execute('''INSERT INTO reference VALUES ('{0}', {1})'''.format(ref, idx))
        self.c.execute("COMMIT")

    def create_indices(self):
        print "Indexing..."
        self.c.execute('''
                       CREATE INDEX align_chrom_position ON align(chrom, position)
                       ''')
