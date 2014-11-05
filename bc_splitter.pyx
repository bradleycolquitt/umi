#! /usr/bin/env python

########################################################################
# Extracts barcode and unique molecular index (UMI) from fastq record. #
# Stores in HDF5 file.
########################################################################
import sys
import os
import gzip
import bmh
import fasta
import pdb
import signal
import sqlite3
import numpy as np
import Levenshtein as lev
from subprocess import Popen

BARCODES = "/home/brad/lib/barcodes/bc2.txt"

class extracter:
    def __init__(self, fastq_fname, dest_fname, compression="blosc"):
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
        chunk_size = 1000000
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

        squal_int = fasta.qual_to_int(qual[self.bc_sub[0]:self.bc_sub[1]])

        # First, check qual. If minimum PHRED is less than 20 anywhere in barcode discard
        if min(squal_int) <= 20: return 0
        # Use Levenshtein distance to identify minimum scoring barcode

        sseq = seq[self.bc_sub[0]:self.bc_sub[1]]
        self.i = 0
        #self.distance_min = 1E6
        #self.val = 0
        #self.ind = 0
        val = 0
        ind = 0
        distance_min = 0
        # for bc in self.barcodes:
        #     self.val = lev.distance(sseq, bc)
        #     if self.val == 0:
        #         self.distance_min = self.val
        #         self.ind = self.i
        #         break
        #     if self.val < self.distance_min:
        #         self.distance_min = self.val
        #         self.ind = self.i
        #     #distance_min = min(ind, lev.distance(sseq, bc))
        #     #self.bc_distances[self.i] = lev.distance(sseq, bc)
        #     self.i += 1
        for bc in self.barcodes:
            val = lev.distance(sseq, bc)
            if val == 0:
                distance_min = val
                ind = self.i
                break
            if val < distance_min:
                distance_min = val
                ind = self.i
            #distance_min = min(ind, lev.distance(sseq, bc))
            #self.bc_distances[self.i] = lev.distance(sseq, bc)
            self.i += 1

        #ind = np.where(self.bc_distances == distance_min)

        #distance_min = self.bc_distances.min()
        #if self.distance_min > 1: return 0
        if distance_min > 1: return 0
        #candidates = np.where(self.bc_distances == distance_min)[0]

        # if len(candidates) > 1:
        #     self.multiple += 1
        #     candidates = self.get_barcode_multiple(candidates, sseq, sqaul_int)
        # if self.distance_min == 1:
        #     #self.max_min += 1
        #     self.ind = self.get_barcode_hard(self.ind, sseq, squal_int)
        # return self.ind
        if distance_min == 1:
            #self.max_min += 1
            ind = self.get_barcode_hard(ind, sseq, squal_int)
        return ind

    # Take candidate barcode from get_barcode
    # Return same barcode index if position(s) contributing to distance have low quality scores
    # Think about this more; may be better just to discard
    # Two scenarios:
    #     1. Poor sequencing causes mismatch: indicated by quality scores
    #     2. Mutation in barcode: not indicated by barcode
    def get_barcode_hard(self, candidate, seq, qual):
        ops = lev.opcodes(seq, self.barcodes[candidate])
        ops = [op for op in ops if op[0]=='replace']
        return candidate

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
        #pdb.set_trace()
        self.c.execute("BEGIN")
        for ind in xrange(len(name)):
            self.c.execute('''INSERT INTO data VALUES ('{0}', '{1}', '{2}', 0, 0)'''.format(name[ind], bc[ind], umi[ind]))
        self.c.execute("COMMIT")
        #self.c.execute("END TRANSACTION")


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
