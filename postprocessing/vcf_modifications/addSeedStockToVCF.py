#!/usr/bin/python
# By: Hans Vasquez-Gross
from __future__ import division
from optparse import OptionParser
import collections
import re
import gzip

parser = OptionParser(usage="usage: %prog [options] output_file\n\nINFO", version="%prog 1.0")
parser.add_option("-m", "--mutantfile", dest="mutant_filename", help="input mutant VCF file name")
parser.add_option("-v", "--verbosity", dest="verbosity_value", type="int", help="Level of verbosity printed to screen; Higher number means more info; 0|1|2; Default: 1", default=1)
(options, args) = parser.parse_args()

if len(args) != 1:
        parser.error("wrong number of arguments")

def addSeedAvailability(lastColumn):

    lastColumnString = "%s;SeedAvailability=%s" % (lastColumn, "Available")

    return lastColumnString

if __name__ == "__main__":
    VERBOSITY = options.verbosity_value
    out_file = args[0]
    out_fh = open(out_file, "w")

    if "gz" in options.mutant_filename:
        in_file_handle = gzip.open(options.mutant_filename, "rb")
    else:
        in_file_handle = open(options.mutant_filename, "r")
    
    for line in in_file_handle:
        line = line.strip()
        lineArr = line.split('\t')
        new_line = line
        if "EFF=" in line:
            seqId = lineArr[0]
            position = lineArr[1]
            lib = lineArr[2]
            ref = lineArr[3]
            mut = lineArr[4]
            lastColumn = lineArr[7]
            newLastColumn = addSeedAvailability(lastColumn)
            new_line = new_line.replace(lastColumn, newLastColumn)
        #print new_line
        out_fh.write(new_line + "\n")


