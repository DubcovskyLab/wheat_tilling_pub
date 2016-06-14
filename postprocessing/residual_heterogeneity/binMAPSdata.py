#!/usr/bin/env python
#Written by Tyson Howell March 11, 2016 at the University of California, Davis

#This script takes a MAPS output file and a contig length file and breaks all chromosomes/scaffolds into bins of a desired size

import argparse
import operator
import math

parser = argparse.ArgumentParser(description='This script takes a MAPS output file and a contig length file and breaks all chromosomes/scaffolds into bins of a desired size')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-l", "--lengths_file", required=True, help="File with length of each chromosome/contig in format <chromName><TAB><chromSize>")
requiredNamed.add_argument("-m", "--maps_file", required=True, help="MAPS file to be binned")
requiredNamed.add_argument("-b", "--bin_size", required=True, help="Bin size to use for splitting in base pairs")
parser.add_argument("-o", "--output", help="name of output file [stdout]")
parser.add_argument("-L", "--output_bin_lengths_file", help="output new lengths file with this name")
args = parser.parse_args()

lengths = open (args.lengths_file)
maps = open (args.maps_file)
bsize = int(args.bin_size)

if args.output:
    out = open(args.output, 'w')
if args.output_bin_lengths_file:
    newlen = open(args.output_bin_lengths_file, 'w')

#Build dictionary of contig lengths
lendict = {}
for line in lengths:
    line = line.rstrip("\r\n")
    fields = line.split("\t")
    lendict[fields[0]] = int(fields[1])
lengths.close()

#Read each line of a MAPS file into a list
mapslist = []
for line in maps:
    if line.startswith("Chrom"): continue #skip header
    line = line.rstrip("\r\n")
    fields = line.split("\t")
    fields[1] = int(fields[1]) #convert position to int for numeric sorting
    mapslist.append(fields)
    
#Sort based first on contig, then on position
mapslist.sort(key = operator.itemgetter(0, 1))
#Loop through list
previous = None
for el in mapslist:
    #determine if length is more than bin size. Skip if not
    if bsize >= lendict[el[0]]:#in this block the contig is smaller than binsize. Just output the line
        if args.output:
            out.write("\t".join(map(str,el)) + "\n")
        else:
            print "\t".join(map(str,el))
        if args.output_bin_lengths_file and previous != el[0]: #If this is the first time we are seeing this contig
            newlen.write(el[0] + "\t" + str(lendict[el[0]]) + "\n")
        previous = el[0]
        continue
    if el[0] != previous: #Check if we are still on the same contig as the last mutation
        #Here we have a new contig with a length greater than bin size.
        #determine number of bins
        numbins = math.ceil(lendict[el[0]]/float(bsize))
        #print el[0], numbins
        
        #Build list of bin sizes
        bins = range(0,int(numbins * bsize) + 1,bsize)
        #print el[0], bins
        
        #remove first bin
        bins = bins[1:-1]
        #print el[0], bins
        
        #Set last bin equal to contig length
        bins.append(lendict[el[0]])
        #print "final", el[0], bins
        curbin = 0
        newbin = 1
        
    #List is sorted, so add MAPS mutations to bins sequentially
    while not el[1] <= bins[curbin]: # if the current mutation is not in the current bin
        curbin += 1 #move to the next bin
        newbin = 1
    if curbin == 0: #For re-naming the contig with bin interval
        prevbin = 1 
    else:
        prevbin = bins[curbin - 1] + 1
    #assign this bin as the "previous" bin before we modify it
    previous = el[0]
    #Update contig name and output
    el[0] = el[0] + "_bin_" + str(prevbin) + "-" + str(bins[curbin])
    if args.output:
        out.write("\t".join(map(str,el)) + "\n")
    else:
        print "\t".join(map(str,el))
    if args.output_bin_lengths_file and newbin:
        newlen.write(el[0] + "\t" + str(bins[curbin] - prevbin + 1) + "\n")
    newbin = 0

    
maps.close()
if args.output: out.close()
if args.output_bin_lengths_file: newlen.close()
