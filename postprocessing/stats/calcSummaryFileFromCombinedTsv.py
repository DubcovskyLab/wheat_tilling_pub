#!/usr/bin/env python
from __future__ import division
import argparse
import pprint
import numpy
import re

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--maps2_file", help="Maps2 results tsv file")
parser.add_argument("-c", "--capture_map", help="Capture map file (in box sync)")
parser.add_argument("-g", "--maps_groupings", help="Maps groupings file (in box sync)")
parser.add_argument("-w", "--wildtype_name", help="Wildtype samplename to skip during averaging;\nSeperate multiple names by ':' (caveat ':' can not be in the name sample name) ", default=None)
parser.add_argument("-i", "--individual_stats", help="Output individual stats into a file with this name")

args = parser.parse_args()


##Main Program
if __name__ == "__main__":
    maps_fh = open(args.maps2_file, "r")

    if(args.wildtype_name):
        wildtype_list = args.wildtype_name.split(':')
    else:
        wildtype_list = list()

    if(args.capture_map):
        captureMapDict = dict()
        captureMap_fh = open(args.capture_map)
        for line in captureMap_fh:
            line = line.strip()
            lineList = line.split('\t')
            captureName = lineList[0]
            sample_name = "Kronos" + lineList[1]
            barcode = lineList[3]
            captureMapDict.setdefault(sample_name, list())
            #captureMapDict[sample_name] = (captureName, barcode)
            captureMapDict[sample_name].append((captureName, barcode))
            #print captureName, sample_name, barcode
    else:
        captureMapDict = dict()

    if(args.maps_groupings):
        mapsGroupingsDict = dict()
        mapsgroupings_fh = open(args.maps_groupings, "r")
        for line in mapsgroupings_fh:
            line = line.strip()
            maps_set = line.split('/')[0]
            filename = line.split('/')[1]
            
            match = re.match("(Kronos\d*)_(KTC\d*_?\d+)?", filename)
            if match:
                sample_name = match.group(1)
                capture = match.group(2)
            else:
                sample_name = ""
                capture = ""
            mapsGroupingsDict.setdefault(sample_name, list())
            #mapsGroupingsDict[sample_name].append((maps_set, capture))
            mapsGroupingsDict[sample_name].append(maps_set)
            
    else:
        mapsGroupingsDict = dict()

    #pprint.pprint(mapsGroupingsDict)

#    ## skip header line
#    header = maps_fh.readline()
#    if "Chrom" not in header:
#        maps_fh.seek(0)
        
    mutantEMSCountDict = dict()
    mutantCountDict = dict()
    mutantHetCountDict = dict()
    mutantWtCovCountDict = dict()
    mutantWtCovHetCountDict = dict()
    mutantWtCovHomCountDict = dict()
    mutantMutCovCountDict = dict()
    mutantMutCovHetCountDict = dict()
    mutantMutCovHomCountDict = dict()
    mutantNonEMStransitionCountDict = dict ()
    mutantNonEMStransitionHetCountDict = dict ()
    mutantNonEMStransitionHomCountDict = dict ()
    for line in maps_fh:
        if line.startswith("Chrom"): continue
        lineList = line.strip().split("\t")
        sample_name = lineList[6]
        bpChange = lineList[10]
        hethom = lineList[7]
        wtcov = int(lineList[8])
        mutcov = int(lineList[9])
        mutantCountDict.setdefault(sample_name, 0)
        mutantCountDict[sample_name] += 1
        mutantWtCovCountDict.setdefault(sample_name, list())
        mutantWtCovCountDict[sample_name].append(wtcov)
        mutantMutCovCountDict.setdefault(sample_name, list())
        mutantMutCovCountDict[sample_name].append(mutcov)
        mutantWtCovHetCountDict.setdefault(sample_name, list())
        mutantWtCovHomCountDict.setdefault(sample_name, list())
        mutantMutCovHetCountDict.setdefault(sample_name, list())
        mutantMutCovHomCountDict.setdefault(sample_name, list())
        if bpChange == "GA" or bpChange == "CT":
        #if bpChange == "GA" or bpChange == "CT" or bpChange == "TC" or bpChange == "AG":
            mutantEMSCountDict.setdefault(sample_name, 0)
            mutantEMSCountDict[sample_name] += 1
        elif bpChange == "AG" or bpChange == "TC":
            mutantNonEMStransitionCountDict.setdefault(sample_name, 0)
            mutantNonEMStransitionCountDict[sample_name] += 1
            if hethom == "het":
                mutantNonEMStransitionHetCountDict.setdefault(sample_name, 0)
                mutantNonEMStransitionHetCountDict[sample_name] += 1
            else:
                mutantNonEMStransitionHomCountDict.setdefault(sample_name, 0)
                mutantNonEMStransitionHomCountDict[sample_name] += 1
        if hethom == "het":
            mutantHetCountDict.setdefault(sample_name, 0)
            mutantHetCountDict[sample_name] += 1
            mutantWtCovHetCountDict[sample_name].append(wtcov)
            mutantMutCovHetCountDict[sample_name].append(mutcov)
        else:
            mutantWtCovHomCountDict[sample_name].append(wtcov)
            mutantMutCovHomCountDict[sample_name].append(mutcov)

    #pprint.pprint(mutantWtCovCountDict)

    if args.individual_stats: ind_fh = open("%s" % args.individual_stats, "w")
    total_wtsnps = 0
    total_mutsnps = 0
    num_wtsamples = 0
    num_mutsamples = 0
    num_wthets = 0
    num_muthets = 0
    num_wtems = 0
    num_mutems = 0
    num_wtnonemstrans = 0
    num_mutnonemstrans = 0
    num_wtnonemstranshet = 0
    num_mutnonemstranshet = 0
    num_wtnonemstranshom = 0
    num_mutnonemstranshom = 0
    #ind_fh.write("num_mutsamples\ttotal_mutsnps\tsample_name\tmutantCountDict[sample_name]\tmutantEMSCountDict[sample_name]\tmutantNonEMStransitionCountDict[sample_name]\tmutantNonEMStransitionHetCountDict[sample_name]\tmutantNonEMStransitionHomCountDict[sample_name]\tmutantHetCountDict[sample_name]\tind_hethom\tems_perc\tavg_wtcov\tavg_mutcov\tavg_wt_to_mutcov\tavg_wthet\tavg_muthet\tavg_wt_to_muthetcov\tavg_wthom\tavg_muthom\tavg_wt_to_muthomcov\n")
    #ind_fh.write("num_mutsamples\ttotal_mutsnps\tsample_name\tmutantCountDict[sample_name]\tmutantEMSCountDict[sample_name]\tmutantNonEMStransitionCountDict[sample_name]\tmutantNonEMStransitionHetCountDict[sample_name]\tmutantNonEMStransitionHomCountDict[sample_name]\tmutantHetCountDict[sample_name]\tind_hethom\tems_perc\tavg_wtcov\tavg_mutcov\tavg_wt_to_mutcov\tavg_wthet\tavg_muthet\tavg_wt_to_muthetcov\tavg_wthom\tavg_muthom\tavg_wt_to_muthomcov\tcapture_set\tbarcode\tmaps_set\n")
    if args.individual_stats: ind_fh.write("num_mutsamples\trunningtotal_mutsnps\tsample_name\ttotal_snps\tems_snps\tnonems_transition\tnonems_transition_hets\tnonems_transition_homs\tcount_hets\tind_hethom\tems_perc\tavg_wtcov\tavg_mutcov\tavg_wt_to_mutcov\tavg_wthet\tavg_muthet\tavg_wt_to_muthetcov\tavg_wthom\tavg_muthom\tavg_wt_to_muthomcov\tcapture_set\tbarcode\tmaps_set\n")
    for sample_name in mutantCountDict:

        if sample_name in wildtype_list:
            mutantCountDict.setdefault(sample_name, 0)
            mutantHetCountDict.setdefault(sample_name, 0)
            mutantEMSCountDict.setdefault(sample_name, 0)
            mutantNonEMStransitionCountDict.setdefault(sample_name, 0)
            if sample_name in mutantCountDict:
                total_wtsnps += mutantCountDict[sample_name]
            if sample_name in mutantHetCountDict:
                num_wthets += mutantHetCountDict[sample_name]
            if sample_name in mutantEMSCountDict:
                num_wtems += mutantEMSCountDict[sample_name]
            if sample_name in mutantNonEMStransitionCountDict:
                num_wtnonemstrans += mutantNonEMStransitionCountDict[sample_name]
            if sample_name in mutantNonEMStransitionHomCountDict:
                num_wtnonemstranshom += mutantNonEMStransitionHomCountDict[sample_name]
            if sample_name in mutantNonEMStransitionHetCountDict:
                num_wtnonemstranshet += mutantNonEMStransitionHetCountDict[sample_name]
            num_wtsamples += 1
        else:
            mutantCountDict.setdefault(sample_name, 0)
            mutantHetCountDict.setdefault(sample_name, 0)
            mutantEMSCountDict.setdefault(sample_name, 0)
            mutantNonEMStransitionCountDict.setdefault(sample_name, 0)
            mutantNonEMStransitionHomCountDict.setdefault(sample_name, 0)
            mutantNonEMStransitionHetCountDict.setdefault(sample_name, 0)
            if sample_name in mutantCountDict:
                total_mutsnps += mutantCountDict[sample_name]
            if sample_name in mutantHetCountDict:
                num_muthets += mutantHetCountDict[sample_name]
            if sample_name in mutantEMSCountDict:
                num_mutems += mutantEMSCountDict[sample_name]
            if sample_name in mutantNonEMStransitionCountDict:
                num_mutnonemstrans += mutantNonEMStransitionCountDict[sample_name]
            if sample_name in mutantNonEMStransitionHomCountDict:
                num_mutnonemstranshom += mutantNonEMStransitionHomCountDict[sample_name]
            if sample_name in mutantNonEMStransitionHetCountDict:
                num_mutnonemstranshet += mutantNonEMStransitionHetCountDict[sample_name]
            num_mutsamples += 1

        if mutantCountDict[sample_name] != 0: ems_perc = mutantEMSCountDict[sample_name] / mutantCountDict[sample_name]
        else: ems_perc = 0
        if mutantHetCountDict[sample_name] and mutantCountDict[sample_name] and mutantHetCountDict[sample_name] and (mutantCountDict[sample_name] - mutantHetCountDict[sample_name]) != 0: ind_hethom = mutantHetCountDict[sample_name] / (mutantCountDict[sample_name] - mutantHetCountDict[sample_name])
        else: ind_hethom = 0
        avg_wtcov = numpy.mean(mutantWtCovCountDict[sample_name])
        avg_mutcov = numpy.mean(mutantMutCovCountDict[sample_name])
        avg_wt_to_mutcov = avg_wtcov / avg_mutcov

        avg_wthetcov = numpy.mean(mutantWtCovHetCountDict[sample_name])
        avg_muthetcov = numpy.mean(mutantMutCovHetCountDict[sample_name])
        avg_wthomcov = numpy.mean(mutantWtCovHomCountDict[sample_name])
        avg_muthomcov = numpy.mean(mutantMutCovHomCountDict[sample_name])
        avg_wt_to_muthetcov = avg_wthetcov / avg_muthetcov
        avg_wt_to_muthomcov = avg_wthomcov / avg_muthomcov

        if sample_name in captureMapDict:
            captureName = "|".join([row[0] for row in captureMapDict[sample_name]])
            barcode = "|".join([row[1] for row in captureMapDict[sample_name]])
        else:
            captureName = ""
            barcode = ""

        if sample_name in mapsGroupingsDict:
            maps_set = ""
            mapsList = mapsGroupingsDict[sample_name]
            if len(mapsList) > 5:
                maps_set = len(mapsList)
            else:
                maps_set = "|".join(mapsList)
            
        else:
            maps_set = ""
        

        if(sample_name not in wildtype_list and args.individual_stats):
            #ind_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (num_mutsamples, total_mutsnps, sample_name, mutantCountDict[sample_name], mutantEMSCountDict[sample_name], mutantNonEMStransitionCountDict[sample_name], mutantNonEMStransitionHetCountDict[sample_name], mutantNonEMStransitionHomCountDict[sample_name], mutantHetCountDict[sample_name], ind_hethom, ems_perc))
            #ind_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (num_mutsamples, total_mutsnps, sample_name, mutantCountDict[sample_name], mutantEMSCountDict[sample_name], mutantNonEMStransitionCountDict[sample_name], mutantNonEMStransitionHetCountDict[sample_name], mutantNonEMStransitionHomCountDict[sample_name], mutantHetCountDict[sample_name], ind_hethom, ems_perc, avg_wtcov, avg_mutcov, avg_wt_to_mutcov, avg_wthetcov, avg_muthetcov, avg_wt_to_muthetcov, avg_wthomcov, avg_muthomcov, avg_wt_to_muthomcov, captureName, barcode))
            ind_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (num_mutsamples, total_mutsnps, sample_name, mutantCountDict[sample_name], mutantEMSCountDict[sample_name], mutantNonEMStransitionCountDict[sample_name], mutantNonEMStransitionHetCountDict[sample_name], mutantNonEMStransitionHomCountDict[sample_name], mutantHetCountDict[sample_name], ind_hethom, ems_perc, avg_wtcov, avg_mutcov, avg_wt_to_mutcov, avg_wthetcov, avg_muthetcov, avg_wt_to_muthetcov, avg_wthomcov, avg_muthomcov, avg_wt_to_muthomcov, captureName, barcode, maps_set))

    if num_wtsamples == 0:
        wtsnps_perline = 0
    else:
        wtsnps_perline = total_wtsnps / num_wtsamples

    mutsnps_perline = total_mutsnps / num_mutsamples
    num_wthoms = total_wtsnps - num_wthets
    num_muthoms = total_mutsnps - num_muthets

    if num_wthoms == 0:
        hethom_wt = 0
    else:
        hethom_wt = num_wthets / num_wthoms

    if num_muthoms == 0:
        hethom_mut = 0
    else:
        hethom_mut = num_muthets / num_muthoms

    if total_wtsnps == 0:
        prcntems_wt = 0
    else:
        prcntems_wt = num_wtems / total_wtsnps

    if total_mutsnps == 0:
        prcntems_mut = 0
    else:
        prcntems_mut = num_mutems / total_mutsnps

    #print "Plant Type(WT/mut)\tTotal # samples\tTotal # SNPs\tAverage # SNPs per line\tTotal # Het SNPs\tHet/Hom ratio\tTotal # non-EMS Transition SNPs\tTotal # EMS SNPs\tPercent EMS SNPs"
    print "Plant Type(WT/mut)\tTotal # samples\tTotal # SNPs\tAverage # SNPs per line\tTotal # Het SNPs\tHet/Hom ratio\tTotal # non-EMS Transition SNPs\tTotal # non-EMS Transition Hets\tTotal # of non-EMS transition Homs\tTotal # EMS SNPs\tPercent EMS SNPs"
    print "wildtype\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (num_wtsamples, total_wtsnps, wtsnps_perline, num_wthets, hethom_wt, num_wtnonemstrans, num_wtnonemstranshet, num_wtnonemstranshom, num_wtems, prcntems_wt)
    print "mutant\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (num_mutsamples, total_mutsnps, mutsnps_perline, num_muthets, hethom_mut, num_mutnonemstrans, num_mutnonemstranshet, num_mutnonemstranshom, num_mutems, prcntems_mut)
    



