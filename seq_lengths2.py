#! /usr/bin/env python

import sys
import gzip
import time

def main(argv):

    file = gzip.open(argv[1])
    total_seq = 0
    lineno = 0

    start = time.time()
    for l in file:
        if lineno % 4 == 1:
            total_seq += len(l)
        lineno += 1

    print "Time elapsed: " + str(time.time() - start)
    print "Total sequence: " + str(total_seq)

if __name__ == "__main__":
    main(sys.argv)
