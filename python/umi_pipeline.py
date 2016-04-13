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
import ConfigParser
from subprocess import Popen
from multiprocessing import Pool

class Meta:
    def __init__(self, meta_fname):
        self.meta_fname = meta_fname
        #self.meta_file = open(meta_fname)
        #self.meta = []

        # Common parameters
        self.data_params = {}
        self.feature_params = {}
        data_params = \
        [
            'data_date',
            'umi_length',
            'min_umi_quality',
            'barcode_set',
            'output',
            'data_files_fname'
        ]
        feature_params = \
        [
            'gtf',
            'group',
            'multimap',
            'multimatch',
            'stranded',
            'min_read_overlap',
        ]
        for p in data_params:
            self.data_params[p] = None

        for p in feature_params:
            self.feature_params[p] = None

        self.data_files = []
        #self.bam_files = []
        # self.params['data_date'] = None
        # self.params['min_umi_quality'] = None
        # self.barcode_set = None
        # self.gtf = None
        # self.allow_multimap = None
        # self.allow_multimatch = None
        # self.stranded = None
        # self.minReadOverlap = None
        # self.output = None
        # self.fastq_files = None
        # self.bam_files = None

        self.parse_config()
        self.load_data_files()
#        self.load_bam_files()

        #pdb.set_trace()
        #header = self.meta_file.readline().strip().split(",")

        ## Parse metadata file
        #i = 0
        # for line in self.meta_file:
        #     if re.search("#", line): continue
        #     sline = line.strip().split(",")
        #     self.meta.append(dict())
        #     j = 0
        #     for element in sline:
        #         self.meta[i][header[j]] = element.strip()
        #         j += 1
        #     i += 1
        #pdb.set_trace()
        #for element in self.meta:

        params_cat = []
        for k,v in self.feature_params.iteritems():
            vs = os.path.basename(os.path.splitext(str(v))[0])
            params_cat.append("-".join([k,vs]))

        params_cat = "_".join(params_cat)
        for data in self.data_files:

            #pdb.set_trace()
            fc_dir = "/".join([os.path.dirname(data['bam']), "featureCounts"])
            if not os.path.exists(fc_dir): os.makedirs(fc_dir)
            #element['bam_prefix'] = os.path.basename(element['bam']).split(".bam")[0]


            data['out_fc'] = "_".join([os.path.basename(data['bam']).split(".bam")[0],
                                       params_cat])

                                          #"overlap" + self.fparams['minReadOverlap']])
            data['out_fc'] = "/".join([fc_dir, data['out_fc']])
            print data['bam'], data['out_fc']

    def parse_config(self):
        config = ConfigParser.ConfigParser()
        config.read(self.meta_fname)

        ## Data-level params
        section = 'data_params'
        options = config.options(section)
        int_options = ['umi_length', 'min_umi_quality']
        string_options = ['data_date', 'barcode_set', 'output', 'data_files_fname']
        bool_options = []

        for option in options:
            if option in int_options:
                self.data_params[option] = config.getint(section, option)
            elif option in string_options:
                self.data_params[option] = config.get(section, option)
            elif option in bool_options:
                self.data_params[option] = config.getboolean(section, option)

        ## Feature-level params
        section = 'feature_params'
        options = config.options('feature_params')
        int_options = ['stranded', 'min_read_overlap']
        string_options = ['gtf', 'group']
        bool_options = ['multimap', 'multimatch']

        for option in options:
            #pdb.set_trace()
            #print option
            if option in int_options:
                self.feature_params[option] = config.getint(section, option)
            elif option in string_options:
                self.feature_params[option] = config.get(section, option)
            elif option in bool_options:
                self.feature_params[option] = config.getboolean(section, option)


    def load_data_files(self):
        data_fname = self.data_params['data_files_fname']
        data_fh = open(data_fname)
        for line in data_fh:
            sline = line.strip().split()
            self.data_files.append({'fastq':sline[0],
                                        'bam':sline[1],
                                    'i7':sline[2]})

    def runFeatureCounts(self):
        print("Running featureCounts...")
        for element in self.data_files:
            print element['bam']
            if os.path.exists(element['out_fc']):
                x = True
                while(True):
                    dec = raw_input("featureCounts file exists. Overwrite [y/n]?")
                    if dec == "n":
                        next
                    elif dec == "y":
                        break
                    else:
                        print "Try again."
            #pdb.set_trace()
            bam_prefix = element['bam'].split(".bam")[0]
            cmd_args = ['featureCounts',
                        '-T', '10',
                        '-R',
                        '-s', str(self.feature_params['stranded']),
                        '-g', self.feature_params['group'],
                        '--minReadOverlap', str(self.feature_params['min_read_overlap']),
                        '-a', self.feature_params['gtf'],
                        '-o', bam_prefix + ".counts"]
            #pdb.set_trace()
            if self.feature_params['multimap']:
                cmd_args = cmd_args + ['-M']
            if self.feature_params['multimatch']:
                cmd_args = cmd_args + ['-O']
            cmd_args = cmd_args + [element['bam']]
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
                    self.data_params['output'],
                    '<',
                    '/home/brad/src/umi/sql/filter.sql']
        cmd_args = " ".join(cmd_args)
#        pdb.set_trace()
        p = Popen(cmd_args, shell=True)
        p.wait()

    def update_info(self):
        import sqlite3
        conn = sqlite3.connect(self.data_params['output'])
        c = conn.cursor()
        fields = self.data_params.keys() + self.feature_params.keys()
        values = self.data_params.values() + self.feature_params.values()
        values = [str(v) for v in values]
        #fields = [" ".join([f, "text"]) for f in fields]

        sql = "DROP TABLE IF EXISTS run_params"
        c.execute(sql)
        sql = "CREATE TABLE run_params (" + ",".join(["`%s` text" for f in fields]) + ")"
        sql = sql % tuple(fields)

        c.execute(sql)

        sql = "INSERT INTO run_params VALUES (" + ",".join(["'%s'" for v in values]) + ")"
        sql = sql % tuple(values)
        c.execute(sql)
        conn.commit()

        sql = "DROP TABLE IF EXISTS input_files"
        c.execute(sql)
        sql = "CREATE TABLE input_files (bam text, fastq text)"
        c.execute(sql)
        sql = "INSERT INTO input_files VALUES (?,?)"
        for element in self.data_files:
            c.execute(sql, (element['bam'], element['fastq']))
        conn.commit()

    def run_bam_hash(self):
        #pdb.set_trace()
        #print("Running bam_hash...")
        for element in self.data_files:
            #pdb.set_trace()
            cmd_args = ['bam_hash2',
                        '-u', str(self.data_params['umi_length']),
                        '-q', str(self.data_params['min_umi_quality']),
                        '--i7', element['i7'],
                        '--barcode', self.data_params['barcode_set'],
                        '--anno', element['out_fc'],
                        '--fastq', element['fastq'],
                        '--outfile', self.data_params['output']]
            cmd_args_join = " ".join(cmd_args)
            try:
                print "Running bam_hash2: " + cmd_args_join
                p = Popen(cmd_args)
                p.wait()
            except OSError as e:
                print e.child_traceback()
                return

def runBamHashWorker(element):
    cmd_args = ['bam_hash2',
                '-u', '10',
                '-q', element['quality'],
                '--i7', element['i7'],
                '--barcode', element['bc'],
                '--anno', element['out_fc'],
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
    for element in obj.data_files:
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
    obj.run_bam_hash()
    #runBamHash(obj)
    obj.update_info()
    obj.collapse()


if __name__ == "__main__":
    main(sys.argv)
0
