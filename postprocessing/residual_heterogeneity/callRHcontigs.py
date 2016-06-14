#!/usr/bin/env python

#This script takes as input a byContig file from detectHeteroZones1.5.py that has been run through addMItoRH.py with the -z option.
#It uses a scoring scheme to determine RH regions, along with some additional criteria
#It outputs two files, one with RH contigs/individuals, and one with non-RH contigs/individuals.

import argparse

parser = argparse.ArgumentParser(description='call RH and non-RH contigs from a file run through detectHeteroZones.py and addMItoRH.py with the -z flag')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-i", "--input", required=True, help="Input file")
requiredNamed.add_argument("-r", "--rh_regions_output", required=True, help="Name of output file for RH regions")
requiredNamed.add_argument("-n", "--non_rh_regions_output", required=True, help="Name of output file for non-RH regions")
parser.add_argument("-t", "--threshold_score",  default=12.5, type=float, help="set score threshold for calling RH [12.5]")
parser.add_argument("-s", "--score_output", action="store_true", help="output an additional column with the calculated score")
args = parser.parse_args()

input = open (args.input)
rhout = open (args.rh_regions_output, "w")
norhout = open (args.non_rh_regions_output, "w")
th = args.threshold_score


for line in input:
    if line.startswith("library"):#reprint header to output files
        rhout.write(line)
        norhout.write(line)
        continue
    score = 0
    line = line.rstrip("\r\n")
    fields = line.split()
    #0-library,1-arm,2-emsCountsPerArmPerLibrary,3-totalmutCountsPerArmPerLibrary,-4prcnt_emsPerArmPerLibrary,-5scaffold_name,
    #6-emsCountsPerContigPerLibrary,7-mutCountsPerContigPerLIbrary,8-contig_prcnt_emsPerLibrary,9-contig_length,10-mutation_densityPerContig,
    #11-psuedomolecule,12-psuedomolecule-start,13-psuedomolecule-end,14-Average-individuals-per-mutation,15-homCountsPerContigPerLibrary
    #Kronos3722      3AL     85      86      0.988372093023  IWGSC_CSS_3AL_scaff_4304353     2       2       1.0     1697    1.17855038303   NA      NA      NA      1.00    2
    muts = int(fields[7])
    density = float(fields[10])
    ems = float(fields[8])
    MI = float(fields[14])
    homperc = float(fields[15])/muts
    #print muts, density, ems, MI, homperc
    
    #Calculate score
    if muts >= 2:                           score += 7
    if ems >= .5 and ems != 1:              score += 2.5
    if ems < .5 and ems != 0:               score += 3.5
    if ems == 0:                            score += 6
    if MI > 1 and MI <= 4:                  score += 0.5
    if MI >4 and MI <6:                     score += 3.5
    if MI >= 6:                             score += 4.5
    if density > 0.2 and density < 0.4:     score += 1.5
    if density >= 0.4 and density < 0.7:    score += 2
    if density >= 0.7:                      score += 3
    #print score
    
    #Check based on homozygosity
    if ems > 0.75 and ems != 1 and homperc >= 0.25 and MI >= 4:  score += th
    if ems <= 0.75 and homperc >= 0.25:                          score += th
    
    #Output to appropriate file
    if score >= th:
        rhout.write(line)
        if args.score_output:
            rhout.write("\t")
            rhout.write(str(score))
        rhout.write("\n")
    if score < th:
        norhout.write(line)
        if args.score_output:
            norhout.write("\t")
            norhout.write(str(score))
        norhout.write("\n")
        
input.close()
rhout.close()
norhout.close()