#!/usr/bin/python
# By: Hans Vasquez-Gross
from __future__ import division
from optparse import OptionParser
import pprint

parser = OptionParser(usage="usage: %prog [options]\n\nPrints to STDOUT the filtered lower confidence file", version="%prog 1.0")
parser.add_option("-v", "--vcfhighconf_file", dest="high_filename", help="input mutant VCF file name with higer confidence")
parser.add_option("-l", "--lowconf_file", dest="low_filename", help="input mutant VCF file name with lower confidence")
(options, args) = parser.parse_args()


if __name__ == "__main__":
    #hc_handle = open(options.high_filename, "r")
    #lc_handle = open(options.low_filename, "r")
    
    ###Skipheader
    #skipline = hc_handle.readline()
    #skipline = hc_handle.readline()
    #skipline = hc_handle.readline()
    #skipline = hc_handle.readline()
    #skipline = hc_handle.readline()
    #skipline = hc_handle.readline()

    #hc_list = list()
    #for line in hc_handle:
    #    line = line.strip()
    #    hc_list.append(line)

    #for line in lc_handle:
    #    line = line.strip()
    #    if line not in hc_list:
    #        print line

    ##NEW
    hc_list = open(options.high_filename).read().strip().split('\n')
    lc_list = open(options.low_filename).read().strip().split('\n')
    lc_list.reverse()
    hc_list.reverse()
    skipline = hc_list.pop()
    skipline = hc_list.pop()
    skipline = hc_list.pop()
    skipline = hc_list.pop()
    skipline = hc_list.pop()
    skipline = hc_list.pop()

    ##Cast list to Set for faster lookup
    hc_set = set(hc_list)
    while lc_list:
        line = lc_list.pop()
        if line not in hc_set:
            print line


