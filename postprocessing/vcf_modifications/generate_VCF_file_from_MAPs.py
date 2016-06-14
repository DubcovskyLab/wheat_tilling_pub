#!/usr/bin/env python

#generate VCF file for Ensembl's VEP tool
#takes mutation file from MAPs (with locus added to it)
#outputs a VCF file for Ensembl's VEP tool

#input python modules
import os, sys
from optparse import OptionParser
from datetime import *

parser = OptionParser(usage="usage: %prog [options] <mut file name> <vcf file name> <vcf indel filename>\n\nINFO", version="%prog 1.0")
parser.add_option("-r", "--referencename", dest="reference_name", help="input reference name for database", default="IWGSC")
parser.add_option("-v", "--verbosity", dest="verbosity_value", type="int", help="Level of verbosity printed to screen; Higher number means more info; 0|1|2; Default: 1", default=1)
(options, args) = parser.parse_args()

if len(sys.argv) != 4:
    raise IOError, 'Usage ./generate_VCF_file_from_MAPs.py <mut file name> <VCF file name> <VCF indel file name>;'

#get filenames
mut_file = open(sys.argv[1])
vcf_file = open(sys.argv[2], 'w')
vcf_indel_file = open(sys.argv[3], 'w')

#get today's date
mydate = date.today()
t = mydate.timetuple()

#write header info for VCF file
vcf_file.write("##fileformat=VCFv4.0\n")
vcf_file.write("##fileDate=" + str(t[0]) + str(t[1]) + str(t[2]) + "\n")
vcf_file.write("##source=Kronos\n")
vcf_file.write("##referece=%s\n" % options.reference_name)
vcf_file.write("##phasing=full\n")
vcf_file.write("##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

#write header infor for VCF indel file
vcf_indel_file.write("##fileformat=VCFv4.0\n")
vcf_indel_file.write("##fileDate=" + str(t[0]) + str(t[1]) + str(t[2]) + "\n")
vcf_indel_file.write("##source=Kronos\n")
#vcf_indel_file.write("##referece=IWGSC\n")
vcf_indel_file.write("##referece=%s\n" % options.reference_name)
vcf_indel_file.write("##phasing=full\n")
vcf_indel_file.write("##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
 

#extract info from mut file and write in vcf file
for line in mut_file:
    if line.startswith("Chrom"): continue
    line_strip = line.strip('\r\n')
    line_split = line_strip.split('\t')
    type_field = line_split[10]
    hethom = line_split[7]
    sample_id = "%s:%s%s%s" % (line_split[14], line_split[2], line_split[1], line_split[5])
    print line_split
    if ("*" in type_field) or ("+" in type_field):
        vcf_indel_file.write(line_split[0] + '\t' + line_split[1] + '\t' + sample_id + '\t' + line_split[2] + '\t' + line_split[5] + '\t40\tPass\tseed_avail=' + line_split[6] + ';hethom=' + hethom + '\n')
    else:
        vcf_file.write(line_split[0] + '\t' + line_split[1] + '\t' + sample_id + '\t' + line_split[2] + '\t' + line_split[5] + '\t40\tPass\tseed_avail=' + line_split[6] + ';hethom=' + hethom + '\n')


print "Done!"
mut_file.close()
vcf_file.close()
vcf_indel_file.close()
