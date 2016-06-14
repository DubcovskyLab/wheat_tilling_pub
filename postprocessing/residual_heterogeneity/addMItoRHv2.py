#!/usr/bin/env python
#Written by Tyson Howell March 8, 2016 at the University of California, Davis

#This script takes as input the byContig output of detectHeteroZonesV1.5.py and the MAPS file used to generate it
#The output is in the byContig file format 

#Version 2: Added option to output number of homozygous mutation calls in an extra column

import argparse
import sys

parser = argparse.ArgumentParser(description='Adds info for average number of individuals each mutation in an RH region occurs to a byContig file')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-r", "--residual_heterogeneity", required=True, help="Specify the '*.byContig' file produced by detectHeteroZonesV1.5.py")
requiredNamed.add_argument("-m", "--maps_file", required=True, help="Specify MAPS file used for generating the RH file")
parser.add_argument("-o", "--output", help="name of output file [stdout]")
parser.add_argument("-z", "--output_hom_counts", action='store_true', help="optionally output another column with the number of heterozygous mutations per individual/contig combination")
#parser.add_argument("-i", "--input_mutation_individual_counts", help="optionally input a file with the number of individuals each mutation occurs in with this name (should be output of -I from a previous run). This can reduce script running time")
args = parser.parse_args()

rh = open (args.residual_heterogeneity)
maps = open (args.maps_file)
if args.output:
    out = open(args.output, 'w')
warnflag = 0

MIdict = {}
mdict = {}
homdict = {}
for line in maps:
    #0-Chrom/Scaffold 1-Pos 2-Ref 3-TotCov 4-WT 5-MA 6-Lib 7-Ho/He 8-WTCov 9-MACov 10-Type 11-LCov 12-#libs 13-InsertType
    #mutation is contig and position
    if line.startswith("Chrom"): continue
    line = line.rstrip("\r\n")
    fields = line.split()
    mutation = fields[0] + " " + fields[1]
    #libtig is library and scaffold
    libtig = fields[6] + " " + fields[0]
    
    #For each contig/position that has a mutation call, build dictionary of how many individuals that mutation occurs in
    if mutation in MIdict:
        MIdict[mutation] += 1
    else:
        MIdict[mutation] = 1    
    
    #Also build dictionary of mutation positions for each individual/contig combination
    if libtig in mdict:
        mdict[libtig].append(fields[1])
    else:
        mdict[libtig] = [fields[1]]
        
    #Also record if mutation is homozygous (1) or heterozygous (0)
    if args.output_hom_counts:
        if fields[7] == 'hom':
            if libtig in homdict:
                homdict[libtig].append(1)
            else:
                homdict[libtig] = [1]
        elif fields[7] == 'het':
            if libtig in homdict:
                homdict[libtig].append(0)
            else:
                homdict[libtig] = [0]
    
maps.close()

#Print header
if args.output_hom_counts:
    if args.output:
        out.write("library	arm	emsCountsPerArmPerLibrary	totalmutCountsPerArmPerLibrary	prcnt_emsPerArmPerLibrary	scaffold_name	emsCountsPerContigPerLibrary	mutCountsPerContigPerLIbrary	contig_prcnt_emsPerLibrary	contig_length	mutation_densityPerContig	psuedomolecule	psuedomolecule-start	psuedomolecule-end	Average-individuals-per-mutation	homCountsPerContigPerLibrary\n")
    else:
        print "library	arm	emsCountsPerArmPerLibrary	totalmutCountsPerArmPerLibrary	prcnt_emsPerArmPerLibrary	scaffold_name	emsCountsPerContigPerLibrary	mutCountsPerContigPerLIbrary	contig_prcnt_emsPerLibrary	contig_length	mutation_densityPerContig	psuedomolecule	psuedomolecule-start	psuedomolecule-end	Average-individuals-per-mutation	homCountsPerContigPerLibrary"
else:
    if args.output:
        out.write("library	arm	emsCountsPerArmPerLibrary	totalmutCountsPerArmPerLibrary	prcnt_emsPerArmPerLibrary	scaffold_name	emsCountsPerContigPerLibrary	mutCountsPerContigPerLIbrary	contig_prcnt_emsPerLibrary	contig_length	mutation_densityPerContig	psuedomolecule	psuedomolecule-start	psuedomolecule-end	Average-individuals-per-mutation\n")
    else:
        print "library	arm	emsCountsPerArmPerLibrary	totalmutCountsPerArmPerLibrary	prcnt_emsPerArmPerLibrary	scaffold_name	emsCountsPerContigPerLibrary	mutCountsPerContigPerLIbrary	contig_prcnt_emsPerLibrary	contig_length	mutation_densityPerContig	psuedomolecule	psuedomolecule-start	psuedomolecule-end	Average-individuals-per-mutation"

#Loop through RH file
for line in rh:
    #0-library	1-arm	2-emsCountsPerArmPerLibrary	3-totalmutCountsPerArmPerLibrary	4-prcnt_emsPerArmPerLibrary	
    #5-scaffold_name	6-emsCountsPerContigPerLibrary	7-mutCountsPerContigPerLIbrary	8-contig_prcnt_emsPerLibrary
    #9-contig_length	10-mutation_densityPerContig	11-psuedomolecule	12-psuedomolecule-start	13-psuedomolecule-end
    if line.startswith("library"): continue
    line = line.rstrip("\r\n")
    fields = line.split()
    libtig = fields[0] + " " + fields[5]
    #For each library/contig file, look up mutation positions in mdict
    if libtig in mdict:
        mutpos = mdict[libtig]
    else:
        sys.stderr.write( "%s is not present in MAPS file!\n" % libtig)
        warnflag = 1
        continue
    #Get number of individuals each mutation occurs in
    inum = 0
    for pos in mutpos:
        #mutation is contig and position
        mutation = fields[5] + " " + pos
        inum += MIdict[mutation]
    #Generate average number of individuals each mutation occurs in for the RH line
    #Calculate and round to two decimal places
    avinum = "{0:.2f}".format(float(inum) / len(mutpos))
    if args.output_hom_counts:
        #print "List: ", len(homdict[libtig]), "Line: ", fields[7]
        homnum = str(sum(homdict[libtig])) #calculate number of homozygous mutation calls
        if args.output:
            out.write(line + "\t" + avinum + "\t" + homnum + "\n")
        else:
            print "%s\t%s\t%s" % (line, avinum, homnum)
    else:
        if args.output:
            out.write(line + "\t" + avinum + "\n")
        else:
            print "%s\t%s" % (line, avinum)
rh.close()

if warnflag == 1:
    sys.stderr.write("WARNING: At least one library-contig set present in the RH file was not present in the MAPS file!\n")

if args.output:
    out.close()