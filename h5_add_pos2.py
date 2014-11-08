#! /usr/bin/env python

###
# Create separate read position H5
###

import os
import pysam
import pdb

import tables as tb


class Alignment(tb.IsDescription):
    name = tb.StringCol(50)
    chrom = tb.UInt8Col()
    position = tb.Int64Col()

class bam_h5:
    def __init__(self, bam_fname, h5_fname):
        self.bam_fname = bam_fname
        self.bam = pysam.Samfile(bam_fname, 'r')
        self.refs = self.bam.references

        self.h5_fname = h5_fname
        self.h5 = 0
        if not os.path.exists(h5_fname):
            self.initialize_h5()
        else:
            self.h5 = tb.open_file(h5_fname, 'a')

    def initialize_h5(self):
        self.h5 = tb.open_file(self.h5_fname, mode = "w")
        filters = tb.Filters(complevel=1, complib="blosc")
        #table = h5.create_table('/', 'data', Alignment, "Alignments", expectedrows=int(size), filters=filters)
        table = self.h5.create_table('/', 'data', Alignment, "Alignments", filters=filters)

    def fill_h5(self):
        table = self.h5.root.data
        record = 0
        for read in self.bam.fetch():
            #pdb.set_trace()
            record = table.row
            record['name'] = read.qname
            record['chrom'] = read.rname
            record['position'] = read.pos
            record.append()
        table.flush()

        self.h5.root.data.cols.name.create_csindex()
        self.h5.root.data.copy("/", "data2", sortby='name', propindexes=True, overwrite=True)
        self.h5.move_node("/", newname="data", name="data2", overwrite=True)

class bc_align:
    def __init__(self, basename):
        self.bc_fname = basename + "_bc.h5"
        self.bc = tb.open_file(self.bc_fname, 'a')
        self.align_fname = basename + "_align.h5"
        self.align = tb.open_file(self.align_fname, 'r')

    def join_bc_align(self):
        bc_data = self.bc.root.data
        align_data = self.align.root.data
        for row in bc_data.itersorted("name"):
            row2 = align_data[row.nrow]
            assert row['name'] == row2['name']
            row['chrom'] = row2['chrom']
            row['position'] = row2['position']
            row.update()
        bc_data.flush()
