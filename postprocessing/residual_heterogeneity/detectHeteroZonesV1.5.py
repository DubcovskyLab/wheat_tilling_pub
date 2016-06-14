#!/usr/bin/python
# By: Hans Vasquez-Gross
# Uses: MAPS2 output TSV file or merged TSV file
# Makes: 3 summary files with prcnt_ems and mutation_density information
# Updated: 2/25/2016
from __future__ import division
from optparse import OptionParser
import numpy
import pprint

parser = OptionParser(usage="usage: %prog [options] output_file_basename\n\nINFO", version="%prog 1.6")
parser.add_option("-m", "--maps2file", dest="maps2_filename", help="input maps2 file")
parser.add_option("-c", "--contigmap", dest="contigmap_filename", help="input contigmap file from Ensembl (OPTIONAL)")
parser.add_option("-t", "--threshold", dest="threshold", type="int", default=1, help="Minimum number of mutations to be considered for the analysis; Default 1")
parser.add_option("-l", "--lengthsfile", dest="lengths_filename", help="input lengths file")
parser.add_option("-w", "--wt", dest="wt_name", help="input wt library identifier")
parser.add_option("-p", "--prcnt_ems", type="float", dest="prcnt", default=1.0)
parser.add_option("-v", "--verbosity", dest="verbosity_value", type="int", help="Level of verbosity printed to screen; Higher number means more info; 0|1|2; Default: 1", default=1)
(options, args) = parser.parse_args()

mutant_types = ["*+",
    "*A",
    "*C",
    "*G",
    "*T",
    "A*",
    "C*",
    "G*",
    "T*",
    "+*",
    "++",
    "+C",
    "+G",
    "+T",
    "+A",
    "A+",
    "C+",
    "G+",
    "T+",
    "AC",
    "AG",
    "AT",
    "CA",
    "CG",
    "GC",
    "GT",
    "TA",
    "TC",
    "TG",
    "CT",
    "GA"
]


def createDataDict(lines):

    dataDictionary = dict()

    for line in lines:
        scaffold = line.split('\t')[0]
        library = line.split('\t')[6]
        hethom = line.split('\t')[7]
        mutantType = line.split('\t')[10]
        if VERBOSITY > 1: print scaffold, library, hethom, mutantType

        dataDictionary.setdefault(library, {key:0 for key in mutant_types})
        dataDictionary.setdefault(options.wt_name, {key:0 for key in mutant_types})
        libraryDict = dataDictionary[library]
        libraryDict[mutantType] += 1
        dataDictionary[library] = libraryDict

    return dataDictionary

def createLengthsDict(in_file):

    lengthsDict = dict()
    in_fh = open(in_file, "r")

    for line in in_fh:
        line = line.strip()
        scaffold = line.split('\t')[0]
        length = line.split('\t')[1]

        lengthsDict[scaffold] = length

    return lengthsDict

def createContigMapDict(in_file):

    contigMapDict = dict()
    in_fh = open(in_file, "r")

    for line in in_fh:
        line = line.strip()
        molecule = line.split('\t')[0]
        start = line.split('\t')[1]
        end = line.split('\t')[2]
        scaffold = line.split('\t')[3]

        contigMapDict[scaffold] = [molecule, start, end]

    return contigMapDict

def generatePercentEMS(out_file_basename):

    dictForLowEMS = dict()
    prcntEMSList = list()

    out_fh = open(out_file_basename + ".byLib.tsv", "w")
    ##Rest of the file
    ##Print to STDOUT just to get a sense of the libraries used
    if VERBOSITY > 0: print "%s\t%s\t%s\t%s\t%s" % ("library", "ems_count", "nonems_count", "total", "prcnt_ems")
    out_fh.write("%s\t%s\t%s\t%s\t%s\n" % ("library", "ems_count", "nonems_count", "total", "prcnt_ems"))
    for library in dataDictionary:
        ems_count = 0
        nonems_count = 0
        total = 0
        for mutantType in dataDictionary[library]:
            if library != options.wt_name:
                if mutantType == 'CT' or mutantType == 'GA':
                    ems_count += dataDictionary[library][mutantType]
                else:
                    nonems_count += dataDictionary[library][mutantType]
                total += dataDictionary[library][mutantType]

        if library != options.wt_name:
            prcnt_ems = ems_count / total
            prcntEMSList.append(prcnt_ems)
            if prcnt_ems <= options.prcnt:
                if VERBOSITY > 0: print "%s\t%s\t%s\t%s\t%s" % (library, ems_count, nonems_count, total, prcnt_ems)
                out_fh.write("%s\t%s\t%s\t%s\t%s\n" % (library, ems_count, nonems_count, total, prcnt_ems))
                dictForLowEMS[library] = (ems_count, nonems_count, total, prcnt_ems)

    median = numpy.median(prcntEMSList)
    mean = numpy.mean(prcntEMSList)
    
    print "Median %EMS for population is ", median
    print "Mean %EMS for population is ", mean
        
    return dictForLowEMS

def findLowEMSArms(out_file_basename):

    mutCountsPerArmDict = dict()
    emsCountsPerArmDict = dict()
    mutCountsPerContigDict = dict()
    emsCountsPerContigDict = dict()
    contigsWithEMS = set()

    ##Process and calculate Mutations and EMS counts per ARM
    for line in lines:
        scaffold = line.split('\t')[0]
        library = line.split('\t')[6]
        hethom = line.split('\t')[7]
        mutantType = line.split('\t')[10]
        
        if scaffold == "3B":
            arm = scaffold
        else:
            arm = scaffold.split("_")[2]

        ##For total mut counts per arm
        mutCountsPerArmDict.setdefault(library, dict())
        mutCountsPerArmDict[library].setdefault(arm, 0)
        mutCountsPerArmDict[library][arm] += 1

        ##For total mut counts per contig
        mutCountsPerContigDict.setdefault(library, dict())
        mutCountsPerContigDict[library].setdefault(arm, dict())
        mutCountsPerContigDict[library][arm].setdefault(scaffold, 0)
        mutCountsPerContigDict[library][arm][scaffold] += 1

        ##For Only EMS muts per arm
        emsCountsPerArmDict.setdefault(library, dict())
        emsCountsPerArmDict[library].setdefault(arm, 0)
        if mutantType == 'CT' or mutantType == 'GA':
            emsCountsPerArmDict[library][arm] += 1

        ##For Only EMS muts per contig
        emsCountsPerContigDict.setdefault(library, dict())
        emsCountsPerContigDict[library].setdefault(arm, dict())
        emsCountsPerContigDict[library][arm].setdefault(scaffold, 0)
        if mutantType == 'CT' or mutantType == 'GA':
            emsCountsPerContigDict[library][arm][scaffold] += 1
            contigsWithEMS.add(scaffold)


    #pprint.pprint(emsCountsPerContigDict)

    ##Write out the file
    out_fh = open(out_file_basename + ".byArm.tsv", "w")
    out_fhContig = open(out_file_basename + ".byContig.tsv", "w")

    out_fh.write("library\tarm\temsCounts\ttotal\tpercent ems\n")
    out_fhContig.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("library", "arm", "emsCountsPerArmPerLibrary", "totalmutCountsPerArmPerLibrary", "prcnt_emsPerArmPerLibrary", "scaffold_name", "emsCountsPerContigPerLibrary", "mutCountsPerContigPerLIbrary", "contig_prcnt_emsPerLibrary","contig_length","mutation_densityPerContig","psuedomolecule","psuedomolecule-start","psuedomolecule-end"))
    savedContigs = set()
    for library in mutCountsPerArmDict:
        for arm in mutCountsPerArmDict[library]:
            prcnt_ems = emsCountsPerArmDict[library][arm] / mutCountsPerArmDict[library][arm]
            if prcnt_ems <= options.prcnt:
                if VERBOSITY > 1: print "%s\t%s\t%s\t%s\t%s" % (library, arm, emsCountsPerArmDict[library][arm], mutCountsPerArmDict[library][arm], prcnt_ems)
                out_fh.write("%s\t%s\t%s\t%s\t%s\n" % (library, arm, emsCountsPerArmDict[library][arm], mutCountsPerArmDict[library][arm], prcnt_ems))
                for scaffold in emsCountsPerContigDict[library][arm]:
                    contig_prcnt_ems = emsCountsPerContigDict[library][arm][scaffold] / mutCountsPerContigDict[library][arm][scaffold]
                    if(mutCountsPerContigDict[library][arm][scaffold] >= options.threshold and contig_prcnt_ems <= options.prcnt):
                        savedContigs.add(scaffold)
                        if VERBOSITY > 1: print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (library, arm, emsCountsPerArmDict[library][arm], mutCountsPerArmDict[library][arm], prcnt_ems, scaffold, emsCountsPerContigDict[library][arm][scaffold], mutCountsPerContigDict[library][arm][scaffold], contig_prcnt_ems)
                        mutationDensity = (int(mutCountsPerContigDict[library][arm][scaffold]) / int(contigLengthsDict[scaffold])) * 1000
                        if options.contigmap_filename and scaffold in contigMapDict:
                            out_fhContig.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (library, arm, emsCountsPerArmDict[library][arm], mutCountsPerArmDict[library][arm], prcnt_ems, scaffold, emsCountsPerContigDict[library][arm][scaffold], mutCountsPerContigDict[library][arm][scaffold], contig_prcnt_ems, contigLengthsDict[scaffold], mutationDensity, contigMapDict[scaffold][0], contigMapDict[scaffold][1], contigMapDict[scaffold][2]))
                        else:
                            out_fhContig.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\n" % (library, arm, emsCountsPerArmDict[library][arm], mutCountsPerArmDict[library][arm], prcnt_ems, scaffold, emsCountsPerContigDict[library][arm][scaffold], mutCountsPerContigDict[library][arm][scaffold], contig_prcnt_ems, contigLengthsDict[scaffold], mutationDensity))


    out_bed = open("lowcontigs.forEnsembl.bed", "w")
    for contig in sorted(savedContigs):
        if "UCW" not in contig:
            out_bed.write("%s\t1\t%s\t%s\n" % (contig, contigLengthsDict[contig], contig))


################
### MAIN PROGRAM
################
if __name__ == "__main__":
    VERBOSITY = options.verbosity_value
    out_file_basename = args[0]
        
    if options.contigmap_filename:
        contigMapDict = createContigMapDict(options.contigmap_filename)

    contigLengthsDict = createLengthsDict(options.lengths_filename)

    lines = [line.rstrip('\n') for line in open(options.maps2_filename)]
    ##Remove header line
    headerline = lines.pop(0)
    if "Chrom" not in headerline:
        lines.insert(0,headerline)

    dataDictionary = createDataDict(lines)
    lowEMSDict = generatePercentEMS(out_file_basename)

    ##Take out lines that do not have lowEMS
    findLowEMSArms(out_file_basename)

