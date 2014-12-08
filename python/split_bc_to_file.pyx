#! /usr/bin/env python

########################################################################
# Extracts barcode and unique molecular index (UMI) from fastq record. #
########################################################################
import sys
import os
import pdb
import gzip
import Levenshtein as lev

import pyximport; pyximport.install()
import bc_util_cy as bc_util
import fasta_cy as fasta

BARCODES = "/home/brad/lib/barcodes/bc8.txt"

## TODO
# filter reads by:
# UMI quality score (min >17)
# discard if >14 A at 3' end
# maximum of 9 G post barcode, discard if otherwise

cdef class extracter:
    cdef public char* fastq_fname
    cdef public char* dest_fname
    cdef int fastq_num_records
    cdef list barcodes
    cdef int bc_end
    cdef list output_files
    cdef tuple bc_sub
    cdef tuple umi_sub
    cdef str sseq

    def __init__(self, fastq_fname):

        ## FASTQ setup ##
        self.fastq_fname = fastq_fname
        print "Buffering fastq data..."

        ## Barcode setup ##
        barcodes_file = open(BARCODES)
        self.barcodes = []
        for l in barcodes_file:
            ls = l.split()
            self.barcodes.append(ls[1])  # expects format of name \t sequence
        self.bc_end = len(self.barcodes)

        ## Output files ##
        if not os.path.exists("split"):
            os.mkdir("split")
        basename = fastq_fname.split(".fastq.gz")[0]
        self.output_files = []

        for idx, item in enumerate(self.barcodes):
            self.output_files.append(open("".join(["split/",basename, "_bc", str(idx+1), ".fastq"]), 'w'))
        self.output_files.append(open("".join(["split/",basename, "_bcUn", ".fastq"]), 'w'))
        ## Other variables ##
        self.bc_sub = (5,11)
        self.umi_sub = (0,5)
        self.sseq = "AAAAAA"

    ## Main function to read in fastq
    cpdef int read_fastq(self):
        cdef int lineno = 0
        cdef int fastqno = 0
        cdef str seq = ""
        cdef str qual = ""
        cdef str name = ""
        cdef int bc = 0
        cdef str umi = ""

        fastq_file = fasta.open_gz(self.fastq_fname)

        print "Reading " + self.fastq_fname
        cdef str l
        for l in fastq_file:
            if lineno % 4 == 0:
                name = l
            if lineno % 4 == 1:
                seq = l

            if lineno % 4 == 3:
                qual = l
                bc = self.get_barcode(seq, qual)
                if bc > -1:
                    self.write_bc(name, seq, qual, bc)
                else:
                    self.write_bc(name, seq, qual, self.bc_end)
            lineno += 1

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
        if (min_qual > 17):
            return seq[sub[0]:sub[1]]
        else:
            return 'N'

    cdef write_bc(self, name, seq, qual, bc):
        self.output_files[bc].write("{0}{1}+\n{2}".format(name, seq, qual))

    cpdef close_files(self):
        [x.close() for x in self.output_files]
