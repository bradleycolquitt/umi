#! /usr/bin/env python

import sys
import argparse

import pyximport; pyximport.install()
import split_bc_to_file as sb

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="fastq")
    args = parser.parse_args()

    ex = sb.extracter(args.fastq)
    ex.read_fastq()
    ex.close_files()

if __name__ == "__main__":
    main(sys.argv)
