#! /usr/bin/env python

###
# Create separate read position H5
###

import os
import pysam
import pdb
import sqlite3
import tables as tb


class Alignment(tb.IsDescription):
    name = tb.StringCol(50)
    chrom = tb.UInt8Col()
    position = tb.Int64Col()

class bam_db:
    def __init__(self, bam_fname, dest_fname):
        self.bam_fname = bam_fname
        self.bam = pysam.Samfile(bam_fname, 'r')
        self.refs = self.bam.references

        self.dest_fname = dest_fname
        self.dest = 0
        self.conn = sqlite3.connect(self.dest_fname)
        self.c = self.conn.cursor()
        self.c.execute('''CREATE TABLE IF NOT EXISTS align (name text, chrom int, position int)''')

    def initialize_dest(self):
        conn = sqlite3.connect(self.dest_fname)
        c = conn.cursor()


        conn.close()

    def fill_dest(self):
        for read in self.bam.fetch():
            self.c.execute('''INSERT INTO align VALUES ('{0}', {1}, {2})'''.format(read.qname, read.rname, read.pos))
        self.conn.commit()


class bc_align:
    def __init__(self, basename):
        self.basename = basename
        self.db = basename + "_bc.db"
        self.conn = sqlite3.connect(self.db)
        self.c = self.conn.cursor()

    def join_bc_align(self):
        self.c.execute('CREATE TABLE IF NOT EXISTS merged (name text, bc int, umi text, chrom, int, position int)')
        self.c.execute('''select data.name,
                          data.bc,
                          data.umi,
                          align.chrom,
                          align.position
                          FROM data, align where align.name = data.name''')
        out = open(self.basename + "_out.txt", 'w')
        for entry in self.c.fetchall():
            out.write("\t".join([str(e) for e in entry]) + "\n")
            #self.c.execute('''INSERT INTO merged VALUES '{0}', {1}, {2}, {3}, {4}'''.format(entry[0].encode(),
             #                                                                         entry[1],
              #                                                                        str(entry[2]),
              #                                                                        entry[3],
              #                                                                        entry[4]))
