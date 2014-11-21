#! /usr/bin/env python

## General pipeline for processing BC/UMI sequences

import os
import pdb
import tables as tb
import bc_splitter2 as bc
import subprocess
import h5_add_pos2 as hap

def main(basename):
    ROOT="/home/brad/src/seq_utils/python/test/"
    BASENAME = basename
    FASTQ=ROOT + BASENAME + ".fastq.gz"
    H5=ROOT + BASENAME
    BAM=FASTQ.split(".")[0]+".bam"
    #COMPS=["blosc", "lzo", "zlib"]
    COMPS=["blosc"]

    ## Align reads to genome using bowtie2
    if not os.path.exists(BAM):
        print "Aligning..."
        cmd_args = "bowtie2 -p 10 -x taeGut1 -U {0} -S {1}; process_sam.sh {1}".format(FASTQ, FASTQ.split(".")[0]+".sam")
        print cmd_args
        p = subprocess.Popen(cmd_args, shell=True)
        p.wait()

    ## Load BC and UMI from fastq into HDF5
    ex = bc.extracter(FASTQ, BASENAME + "_bc.h5")
    ret = ex.read_fastq()
    ex.h5.close()

    ## Update HDF5 with read aligned chrom and pos using read name as key
    hap1 = hap.bam_h5(BAM, BASENAME + "_align.h5")
    print "fill_h5"
    hap1.fill_h5()
    hap1.h5.close()

    print "join"
    ba = hap.bc_align(BASENAME)
    ba.join_bc_align()
