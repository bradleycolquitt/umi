#! /usr/bin/env python

import os
import pysam
import pdb
import sqlite3

cdef class bam_db:
    cdef str bam_fname
    cdef object bam
    cdef tuple refs
    cdef long bam_counts
    cdef str dest_fname
    cdef object conn
    cdef object c

    def __init__(self, bam_fname, dest_fname):
        self.bam_fname = bam_fname
        self.bam = pysam.Samfile(bam_fname, 'r')
        self.refs = self.bam.references
        self.bam_counts = self.bam.count()

        self.dest_fname = dest_fname

        self.conn = sqlite3.connect(self.dest_fname, isolation_level=None)
        self.c = self.conn.cursor()

    def fill_join(self):
        self.fill_db()
        self.join_bc_align()

    cpdef fill_db(self):
        # Insert BAM data into db
        print "Filling " + self.dest_fname
        self.c.execute('''CREATE TABLE IF NOT EXISTS align (name text,
                                                            chrom int,
                                                            position int,
                                                            strand int)''')
        cdef long readno = 0
        chunk_size = self.bam_counts

        cdef object read
        cdef int pos
        self.c.execute("BEGIN TRANSACTION")
        for read in self.bam.fetch():
            if not read.is_unmapped:
                if read.is_reverse:
                    pos = read.alen
                else:
                    pos = read.pos
                self.c.execute('''INSERT INTO align VALUES (?, ?, ?, ?)''', (read.qname,
                                                                             read.rname,
                                                                             pos,
                                                                             read.is_reverse))
                readno += 1
            if readno == chunk_size:
                self.c.execute("COMMIT")
                self.c.execute("BEGIN TRANSACTION")
                readno = 0
        self.c.execute("COMMIT")

        # Create index for name, chrom, position
        print "Indexing..."
        self.create_indices()

        # Create references table containing key for refid and human-readable chromosome label
        print "Create Ref id table"
        self.c.execute("BEGIN TRANSACTION")
        self.c.execute('''CREATE TABLE IF NOT EXISTS reference (name text, chrom int)''')
        for idx, ref in enumerate(self.refs):
            self.c.execute('''INSERT INTO reference VALUES ('{0}', {1})'''.format(ref, idx))
        self.c.execute("COMMIT")

    cpdef join_bc_align(self):
        self.c.execute('''select data.name,
                          data.bc,
                          data.umi,
                          align.chrom,
                          align.position
                          FROM data INNER JOIN align where align.name = data.name''')

        out = open(self.dest_fname + "_out.txt", 'w')
        cdef tuple entry
        for entry in self.c.fetchall():
            out.write("\t".join([str(e) for e in entry]) + "\n")

    def create_indices(self):
        self.c.execute('''
                       CREATE INDEX align_name ON align(name)
                       ''')
        self.c.execute('''
                       CREATE INDEX align_chrom_position ON align(chrom, position)
                       ''')
