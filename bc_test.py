#! /usr/bin/env python

import bc_splitter1 as bc

def main():
    fastq_fname = "/home/brad/data2/fastq/141016/Lane1/index_1_fastq/L1_AAGGGA_L001_R1_001.fastq.gz"
    h5_fname = "/home/brad/src/seq_utils/python/test2.h5"
    ex = bc.extracter(fastq_fname, h5_fname)
    ex.read_fastq()



if __name__ == "__main__":
    main()
