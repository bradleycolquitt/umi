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
import numpy as np
import Levenshtein as lev
from tables import *

BARCODES = "/home/brad/lib/barcodes/bc2.txt"

class extracter:
    def __init__(self, fastq_fname, h5_fname, compression="blosc"):
        try:
            signal.signal(signal.SIGINT, self.signal_term_handler)
            ## FASTQ setup ##
            self.fastq_fname = fastq_fname
            print "Buffering fastq data..."

            self.fastq_file = fasta.open_gz(fastq_fname)

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

            ## HDF5 setup ##
            self.h5_fname = h5_fname
            self.h5 = 0

            # initialize h5 if does not exist
            if not os.path.exists(self.h5_fname):
                self.h5 = initialize_h5(self.h5_fname, self.fastq_num_records, compression)
            else:
                self.h5 = open_file(h5_fname, 'a')
            self.table = self.h5.root.data

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
            self.bc_distances = np.zeros(self.num_barcodes)
        except:
            self.h5.close()

    def read_fastq(self):
        lineno = 0
        name = ""
        seq = ""
        qual = ""
        bc = ""
        umi = ""

        print "Reading " + self.fastq_fname
        for l in self.fastq_file:
            if lineno % 4 == 0:
                name = l.split()[0].strip('@')
            if lineno % 4 == 1:
                seq = l
            if lineno % 4 == 3:
                qual = l
                bc = self.get_barcode(seq, qual)
                umi = self.get_umi(seq, qual)
                self.insert_to_h5(name, bc, umi)
            lineno += 1
        self.h5.root.data.flush()

        print "Sort..."
        self.h5.root.data.cols.name.create_csindex()
        self.h5.root.data.copy("/", "data2", sortby='name', propindexes=True)
        self.h5.move_node("/", newname="data", name="data2", overwrite=True)
        #self.h5.close()

    # Search sequence [5:11] for barcodes
    # Returns index of barcode with minimum distance
    def get_barcode(self, seq, qual):

        squal_int = fasta.qual_to_int(qual[self.bc_sub[0]:self.bc_sub[1]])

        # First, check qual. If minimum PHRED is less than 20 anywhere in barcode discard
        if min(squal_int) <= 20: return 0
        # Use Levenshtein distance to identify minimum scoring barcode

        sseq = seq[self.bc_sub[0]:self.bc_sub[1]]
        self.i = 0

        for bc in self.barcodes:
            self.bc_distances[self.i] = lev.distance(sseq, bc)
            self.i += 1

        distance_min = self.bc_distances.min()
        if distance_min > 1: return 0
        candidates = np.where(self.bc_distances == distance_min)[0]

        # if len(candidates) > 1:
        #     self.multiple += 1
        #     candidates = self.get_barcode_multiple(candidates, sseq, sqaul_int)
        if distance_min == 1:
            #self.max_min += 1
            candidates = self.get_barcode_hard(candidates, sseq, squal_int)
        return candidates

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

    def insert_to_h5(self, name, bc, umi):
        table = self.table
        record = table.row
        record['name'] = name
        record['bc'] = bc
        record['umi'] = umi
        record.append()

    def signal_term_handler(self, signal, frame):
        self.h5.close()

class Read(IsDescription):
    name = StringCol(50)
    bc = UInt8Col()
    umi = StringCol(5)
    chrom = UInt8Col()
    position = Int64Col()

# include functionality to read aligner info from MySQL db and add as attributes
#    * aligner
#    * alignment date
#    * command line
#    * barcode file
def initialize_h5(h5_fname, size, compression):
    h5 = open_file(h5_fname, mode = "w")
    filters = Filters(complevel=1, complib=compression)
    table = h5.create_table('/', 'data', Read, "Reads", expectedrows=int(size), filters=filters)

    # Set attributes
    table.attrs.bc_file = BARCODES

    return h5
