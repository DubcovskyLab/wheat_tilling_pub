#!/usr/bin/env python
from __future__ import division
from itertools import izip
import argparse
import pprint

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--assay_file", help="assay file from Mapspart1 script")
parser.add_argument("-c", "--coverage_minimum", type=int, help="minimum coverage for a baseposition (hetmincov value)")
parser.add_argument("-w", "--wildtype_name", help="Wildtype samplename to skip during averaging;\nSeperate multiple names by ':' (caveat ':' can not be in the name sample name) ", default="")
args = parser.parse_args()

##Function Definitions
def initDictFromHeader(header):

    dictionary = dict()
    sampleOrderList = list()
    headerList = header.split("\t")
    #print headerList
    for i in range(1, (len(headerList) - 1)):
        sample_name = headerList[i]
        #print headerList[i]
        dictionary.setdefault(sample_name, 0)
        sampleOrderList.append(sample_name)

    return (dictionary, sampleOrderList)

##Main Program
if __name__ == "__main__":
    if(args.assay_file):
        assay_fh = open(args.assay_file, "r")
    else:
        raise Exception('Assay file opened failed!')

    if(args.coverage_minimum):
        minCov = int(args.coverage_minimum)
    else:
        raise Exception('Coverage minimum not set')

    if(args.wildtype_name):
        wildtype_list = args.wildtype_name.split(':')
    else:
        wildtype_list = list()

    header = assay_fh.readline()
    #print header

    ##Parse coverage for a given minCov value
    (covDict, sampleOrderList) = initDictFromHeader(header.strip()) 
    #print 'Init dict'
    #pprint.pprint(covDict)
    current_coverage = 0
    for line in assay_fh:
        current_coverage += 1
        lineList = line.strip().split("\t")
        for my_tuple in zip(sampleOrderList, range(1, (len(lineList) - 1))):
            sample_name = my_tuple[0]
            i = my_tuple[1]
            coverage = int(round(float(lineList[i])))
            if current_coverage >= minCov:
                covDict[sample_name] += coverage

    assay_fh.close() 

    #print 'data covDict'
    #pprint.pprint(covDict)
    print "Sample_name\tsequence_space_at_mincov"
    mutCovSum = 0
    count = 0
    for sample_name in sorted(covDict):
        count += 1
        print "%s\t%s\t%s" % (count, sample_name, covDict[sample_name])

        ##Skip user specified wildtype name
        if sample_name not in wildtype_list:
            mutCovSum += covDict[sample_name]

    num_samples = ""
    if args.wildtype_name:
        num_skipped_samples = len(wildtype_list)
        mutCovAvg = mutCovSum / (len(covDict) - num_skipped_samples)
        num_samples = len(covDict) - num_skipped_samples
    else:
        mutCovAvg = mutCovSum / len(covDict)
        num_samples = len(covDict)
    print "%s\t%s\t%.9f" % ("Average",  num_samples, mutCovAvg)
    

