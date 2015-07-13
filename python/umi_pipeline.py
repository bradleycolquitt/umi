#! /usr/bin/env python

"""
read in metadata file containing:
    date
    i7
    read1 fastq.gz
    read2 bam
    gtf
    output

featureCounts
bam_hash2
"""

import sys
import re
import argparse
import pdb
from subprocess import Popen
from multiprocessing import Pool

class Meta:
    def __init__(self, meta_fname):
        self.meta_file = open(meta_fname)
        self.meta = []
        header = self.meta_file.readline().strip().split(",")


        i = 0
        for line in self.meta_file:
            if re.search("#", line): continue
            sline = line.strip().split(",")
            self.meta.append(dict())
            j = 0
            for element in sline:
                self.meta[i][header[j]] = element
                j += 1
            i += 1

    def runFeatureCounts(self):
        for element in self.meta:

            bam_prefix = element['bam'].split(".bam")[0]
            cmd_args = ['featureCounts',
                        '-T', '10',
                        '-R',
                        '-s', '1',
                        '-a', element['gtf'],
                        '-o', bam_prefix + ".counts",
                        element['bam']]
            cmd_args_join = " ".join(cmd_args)
            try:
                print "Running featureCounts: " + cmd_args_join
                p = Popen(cmd_args)
                p.wait()
            except:

                return

def runBamHashWorker(element):

    cmd_args = ['bam_hash2',
                '-q', '20',
                '--i7', element['i7'],
                '--barcode', element['bc'],
                '--anno', element['bam']+'.featureCounts',
                '--bam', element['bam'],
                '--fastq', element['fastq'],
                '--outfile', element['output']]
    cmd_args_join = " ".join(cmd_args)
    try:
        print "Running bam_hash2: " + cmd_args_join
        p = Popen(cmd_args)
        p.wait()
    except OSError as e:
        print e.child_traceback()
        return

def runBamHash(obj):
    pool = Pool(processes=10)
    # pdb.set_trace()
    for element in obj.meta:
        pool.apply_async(runBamHashWorker, (element, ))
    pool.close()
    pool.join()


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', dest='meta')
    args = parser.parse_args()

    obj = Meta(args.meta)
    #pdb.set_trace()
#    obj.runFeatureCounts()
    runBamHash(obj)


if __name__ == "__main__":
    main(sys.argv)
