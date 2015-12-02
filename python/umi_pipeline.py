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


import os
import sys
import re
import shutil
import argparse
import pdb
from subprocess import Popen
from multiprocessing import Pool

class Meta:
    def __init__(self, meta_fname):
        self.meta_file = open(meta_fname)
        self.meta = []
        header = self.meta_file.readline().strip().split(",")

        ## Parse metadata file
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

        for element in self.meta:
            #pdb.set_trace()
            fc_dir = "/".join([os.path.dirname(element['bam']), "featureCounts"])
            if not os.path.exists(fc_dir): os.makedirs(fc_dir)
            #element['bam_prefix'] = os.path.basename(element['bam']).split(".bam")[0]
            element['out_fc'] = "-".join([os.path.basename(element['bam']).split(".bam")[0],
                                          os.path.basename(element['gtf']).split(".gtf")[0],
                                          "overlap" + element['minReadOverlap']])
            element['out_fc'] = "/".join([fc_dir, element['out_fc']])

    def runFeatureCounts(self):
        for element in self.meta:
            if os.path.exists(element['out_fc']):
                while(True):
                    dec = raw_input("featureCounts file exists. Overwrite [y/n]?")
                    if dec == "n":
                        return
                    elif dec == "y":
                        break
                    else:
                        print "Try again."

            bam_prefix = element['bam'].split(".bam")[0]
            cmd_args = ['featureCounts',
                        '-T', '10',
                        '-R',
#                        '-t', 'gene',
                        '-s', '1',
                        '-g', 'gene_id',
                        '-M',
#                        '-O',
                        '--minReadOverlap', element['minReadOverlap'],
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

            shutil.move(element['bam'] + ".featureCounts",
                        element['out_fc'])
            shutil.move(element['bam'].split(".bam")[0] + ".counts.summary",
                        element['out_fc'] + ".summary")

def runBamHashWorker(element):
    cmd_args = ['bam_hash2',
                '-u', '10',
                '-q', element['quality'],
                '--i7', element['i7'],
                '--barcode', element['bc'],
                '--anno', element['out_fc'],
                #'--bam', element['bam'],
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
    for element in obj.meta:
        pool.apply(runBamHashWorker, (element, ))
       #pool.apply_async(runBamHashWorker, (element, ))
    pool.close()
    pool.join()

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', dest='meta')
    args = parser.parse_args()

    obj = Meta(args.meta)
    obj.runFeatureCounts()
    runBamHash(obj)


if __name__ == "__main__":
    main(sys.argv)
