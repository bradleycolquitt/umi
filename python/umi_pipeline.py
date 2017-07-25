#! /usr/bin/env python

import sys
import argparse
import umi_util as uu

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('process', help="Process to run. 'pipeline', 'split'")
    parser.add_argument('meta', help="TSV file containing samples to process")
    args = parser.parse_args()

    obj = uu.Meta(args.meta)

    if args.process=="pipeline":
        obj.runFeatureCounts()
        obj.run_bam_hash()
        obj.update_info()
        obj.collapse()
    elif args.process=="split":
        obj.run_bam_split()

if __name__ == "__main__":
    main(sys.argv)
