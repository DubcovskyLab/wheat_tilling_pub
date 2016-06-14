#!/usr/bin/python
# By: Hans Vasquez-Gross
from __future__ import division
from optparse import OptionParser
import collections
import re
import pprint

parser = OptionParser(usage="usage: %prog [options] output_file\n\nINFO", version="%prog 1.0")
parser.add_option("-m", "--mutantfile", dest="mutant_filename", help="input mutant VCF file name")
parser.add_option("-r", "--residheterofile", dest="hetero_filename", help="input residual heterogeneity file name")
parser.add_option("-s", "--skip_samples", dest="skip_samples", help="input sample names to skip when merging (seperated by ':' characters", default="")
parser.add_option("-v", "--verbosity", dest="verbosity_value", type="int", help="Level of verbosity printed to screen; Higher number means more info; 0|1|2; Default: 1", default=1)
(options, args) = parser.parse_args()

#if len(args) != 1:
#        parser.error("wrong number of arguments")

if __name__ == "__main__":
    skip_samples_list = options.skip_samples.split(':')

    ##Proces RH file
    hetero_file_handle = open(options.hetero_filename, "r")
    heteroDataDict = dict()
    for line in hetero_file_handle:
        line = line.strip()
        lineArr = line.split('\t')
        library = lineArr[0]
        seqId = lineArr[5]
        ##IWGSC_CSS_3AL_scaff_4330770_bin_1-10000
        if "bin" in seqId:
            temp = seqId.split('bin_')[1]
            startPos = temp.split('-')[0]
            endPos = temp.split('-')[1]
            seqId = seqId.split('_bin')[0]
            rangeLoc = (startPos, endPos)
        else:
            rangeLoc = (0,0)
        #heteroDataDict.setdefault(seqId, list())
        #heteroDataDict[seqId].append(library)
        heteroDataDict.setdefault(seqId, dict())
        heteroDataDict[seqId].setdefault(library, list())
        heteroDataDict[seqId][library].append(rangeLoc)

    #print "Printing heteroDataDict"
    #pprint.pprint(heteroDataDict)

    in_file_handle = open(options.mutant_filename, "r")

    ##Skipheader
    skipline = in_file_handle.readline().strip()
    print skipline
    skipline = in_file_handle.readline().strip()
    print skipline
    skipline = in_file_handle.readline().strip()
    print skipline
    skipline = in_file_handle.readline().strip()
    print skipline
    skipline = in_file_handle.readline().strip()
    print skipline
    skipline = in_file_handle.readline().strip()
    print skipline


    for line in in_file_handle:
        line = line.strip()
        lineArr = line.split('\t')
        seqId = lineArr[0]
        pos = lineArr[1]
        library = lineArr[2].split(':')[0]
        info = lineArr[7]

        if library not in skip_samples_list:
            if seqId in heteroDataDict:
                if library in heteroDataDict[seqId]:
                    for posTuple in heteroDataDict[seqId][library]:
                        ##IF the tuple has a defined region, only tag ones that occur within the range, if not print the normal maps line
                        if posTuple[0] > 0:
                            if pos >= posTuple[0] and pos <= posTuple:
                                print line + ";resid_hetero=\"<b><font color='red'>WARNING:</font></b> This SNP is in a region of residual heterogeneity and may be a natural variant."
                            else:
                                print line
                        else:
                            print line + ";resid_hetero=\"<b><font color='red'>WARNING:</font></b> This SNP is in a region of residual heterogeneity and may be a natural variant."
                else:
                    print line
            else:
                print line
        else:
            print line

