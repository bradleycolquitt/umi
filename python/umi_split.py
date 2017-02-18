#! /usr/bin/env python

import sys
import argparse
import umi_util as uu

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('meta', help="CSV file containing samples to process")
    args = parser.parse_args()

    obj = uu.Meta(args.meta)
    obj.run_bam_split()

if __name__ == "__main__":
    main(sys.argv)
