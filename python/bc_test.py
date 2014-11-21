#! /usr/bin/env python

import unittest

import sys
import sqlite3
import traceback as tb
from subprocess import Popen

import pyximport; pyximport.install()
import bc_splitter_cy as bc
import bam_sql_cy as bs

path = "/home/brad/src/seq_utils/test_files/"
base = "test_dna_40000"
test_fastq = "".join([path, base, ".fastq.gz"])
test_bam = "".join([path, base, ".bam"])
test_db = "".join([path, base, ".db"])

class BarcodeSplitterTest(unittest.TestCase):
    def setUp(self):
        try:
            self.ex=bc.extracter(test_fastq, test_db)
            self.bs1 = bs.bam_db(test_bam, test_db)
        except:
            self.fail( tb.print_exc())

    def test_read_fastq_sql(self):
        try:
            self.ex.read_fastq_sql()
        except:
            self.fail(tb.print_exc())

    def test_bam_sql(self):
        try:
            self.bs1.fill_db()
            conn = sqlite3.connect(test_db)
            c = conn.cursor()
            c.execute("SELECT * FROM align LIMIT 0, 10")
            for res in c.fetchall():
                print res
        except:
            self.fail(tb.print_exc())
        finally:
            print "Removing test db..."
            cmd_args = ["rm", test_db]
            p = Popen(cmd_args)
            p.wait()


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(BarcodeSplitterTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
