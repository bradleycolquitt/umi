#! /usr/bin/env python

import numpy as np

def qual_to_int(qual):
    return np.array([ord(c) - 33 for c in qual])
