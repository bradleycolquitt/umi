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
        #pdb.set_trace()
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
                self.meta[i][header[j]] = element.strip()
                j += 1
            i += 1
        #pdb.set_trace()
        for element in self.meta:
            #pdb.set_trace()
            fc_dir = "/".join([os.path.dirname(element['bam']), "featureCounts"])
            if not os.path.exists(fc_dir): os.makedirs(fc_dir)
            #element['bam_prefix'] = os.path.basename(element['bam']).split(".bam")[0]
            element['out_fc'] = "-".join([os.path.basename(element['bam']).split(".bam")[0],
                                          os.path.basename(element['gtf']).split(".gtf")[0],
                                          "overlap" + element['minReadOverlap']])
            element['out_fc'] = "/".join([fc_dir, element['out_fc']])
            print element['bam'], element['out_fc']

    def read_config(self):
        config = ConfigParser.ConfigParser()
        config.read('config2.config')

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
            #pdb.set_trace()
            bam_prefix = element['bam'].split(".bam")[0]
            cmd_args = ['featureCounts',
                        '-T', '10',
                        '-R',
#                        '-t', 'gene',
                        '-s', '1',
                        '-g', 'gene_id',
#                        '-M',
#                        '-O',
                        '--minReadOverlap', element['minReadOverlap'],
                        '-a', element['gtf'],
                        '-o', bam_prefix + ".counts",
                        element['bam']]
            cmd_args_join = " ".join(cmd_args)
            try:
                print "Running featureCounts: " + cmd_args_join
                errlog = open(element['out_fc'] + ".errlog", 'w')
                p = Popen(cmd_args, stderr=errlog)
                p.wait()
                errlog.close()
            except:
                return

            shutil.move(element['bam'] + ".featureCounts",
                        element['out_fc'])
            shutil.move(element['bam'].split(".bam")[0] + ".counts.summary",
                        element['out_fc'] + ".summary")

    def collapse(self):
        cmd_args = ['sqlite3',
                    self.meta[0]['output'],
                    '<',
                    '/home/brad/src/umi/python/filter.sql']
        p = Popen(cmd_args, shell=T)
        p.wait()

    def update_info(self):
        import sqlite3
        conn = sqlite3.connect(self.meta[0]['output'])
        c = conn.cursor()
        fields = self.meta[0].keys()
        fields = [" ".join([f, "text"]) for f in fields]

        sql = "DROP TABLE IF EXISTS analysis_info"
        c.execute(sql)
        sql = "CREATE TABLE analysis_info (" + ",".join(fields) + ")"
        c.execute(sql)
        sql = "INSERT INTO analysis_info VALUES (" + ",".join(["'?'" for f in fields]) + ")"
        for element in self.meta:
            #sql1 = sql tuple(element.values())
            c.execute(sql, tuple(element.values()) )
        conn.commit()

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
    parser.add_argument('meta', help="CSV file containing samples to process")
    args = parser.parse_args()

    obj = Meta(args.meta)
    obj.runFeatureCounts()
    runBamHash(obj)
    obj.update_info()
    obj.collapse()


if __name__ == "__main__":
    main(sys.argv)
