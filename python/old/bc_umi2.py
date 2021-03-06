#! /usr/bin/env python

## General pipeline for processing BC/UMI sequences

import os
import tables as tb
import bc_splitter2 as bc
import subprocess
import h5_add_pos as hap

ROOT="/home/brad/src/seq_utils/python/test/"
FASTQ=ROOT + "test_dna3.fastq.gz"
H5=ROOT + "test_dna3.h5"
BAM=FASTQ.split(".")[0]+".bam"

def main():

    ## Load BC and UMI from fastq into HDF5
    if not os.path.exists(H5):
        ex = bc.extracter(FASTQ, H5)
        ret = ex.read_fastq()
        ex.h5.close()

    ## Align reads to genome using bowtie2
    if not os.path.exists(BAM):
        print "Aligning..."
        cmd_args = "bowtie2 -p 10 -x taeGut1 -U {0} -S {1}; process_sam.sh {1}".format(FASTQ, FASTQ.split(".")[0]+".sam")
        print cmd_args
        p = subprocess.Popen(cmd_args, shell=True)
        p.wait()

    ## Update HDF5 with read aligned chrom and pos using read name as key
    hap1 = hap.bam_h5(BAM, H5)
    print "Indexing read column, create_csindex()..."
    hap1.h5.root.data.cols.name.create_csindex()
    print "update()..."
    hap1.update()

    hap1.h5.root.data.cols.name.remove_index()
    print "Indexing read column, create_index()..."
    hap1.h5.root.data.cols.name.create_index()
    print "update()..."
    hap1.update()

    hap1.h5.root.data.cols.name.remove_index()
    print "Indexing read column, create_index(optlevel=1, kind='ultralight')..."
    hap1.h5.root.data.cols.name.create_index(optlevel=1, kind='ultralight')
    print "update()..."
    hap1.update()
