#! /usr/bin/env python

import numpy as np
import subprocess

import cStringIO
io_method = cStringIO.StringIO

def qual_to_int(qual):
    return np.array([ord(c) - 33 for c in qual])

def open_gz(fname):
    p = subprocess.Popen(["zcat", fname], stdout = subprocess.PIPE)
    return io_method(p.communicate()[0])

def count_records(fname):
    file = open_gz(fname)
    count = 0
    for l in file:
        count += 1
    return count / 4
