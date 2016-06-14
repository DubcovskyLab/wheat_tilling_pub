#!/usr/bin/python
# By: Hans Vasquez-Gross
from __future__ import division
from optparse import OptionParser
import collections
import re
import pprint

parser = OptionParser(usage="usage: %prog [options]\n\nThis program prints a summary to STDOUT and creates a stat file named individuals_stats.tsv", version="%prog 1.0")
parser.add_option("-m", "--mutantfile", dest="mutant_filename", help="input mutant VEP VCF file name")
parser.add_option("-s", "--skip_samples", dest="skip_samples", help="input sample names to skip when merging (seperated by ':' characters", default="")
parser.add_option("-v", "--verbosity", dest="verbosity_value", type="int", help="Level of verbosity printed to screen; Higher number means more info; 0|1|2; Default: 1", default=1)
(options, args) = parser.parse_args()

#if len(args) != 1:
#        parser.error("wrong number of arguments")

if __name__ == "__main__":
    skip_samples_list = options.skip_samples.split(':')

    if options.mutant_filename:
        in_file_handle = open(options.mutant_filename, "r")
    else:
        raise Exception('VEP VCF mutant file not provided')
        


##dataDict dictionary is a key of consequences with the value being a list of gene models
    dataDict = dict()
    for line in in_file_handle:

        ##Skip header
        if line.startswith("#"):
            continue

        line = line.strip()
        lineArr = line.split('\t')
        seqId = lineArr[0]
        pos = lineArr[1]
        library = lineArr[2].split(':')[0]
        info = lineArr[7]

        if library not in skip_samples_list:
            infoList = info.split(';')
            csq = ""
            for item in infoList:
                if "CSQ" in item or "VE" in item:
                    csq = item

            if "CSQ" in csq or "VE" in csq:
                for csq in csq.split(','):
                    csqArr = csq.split('|')
                    model = csqArr[1]
                    transcript = csqArr[2]
                    seqtype = csqArr[3]
                    conseq_type = csqArr[4]
                    for conseq in conseq_type.split('&'):
                        dataDict.setdefault(conseq, list())
                        dataDict[conseq].append(model)


    #pprint.pprint(dataDict)

    total = 0
    gm_list = list()
    for conseq in sorted(dataDict):
        count = len(dataDict[conseq])
        gm_list.extend(dataDict[conseq])
        print "%s\t%s" % (conseq, count)
        total += count

    print "%s\t%s" % ("total", total)

    ##Get uniq list of gene models
    ##then Filter out empty strings
    uniq_list = list(set(gm_list))
    filtered_list = filter(None, uniq_list)
    gm_list = filtered_list
    print "\nTotal Gene Models with mutations\t%s" % len(filtered_list) 

    ##find the number of gene models with premature stop codons
    premature_stop_list = dataDict["stop_gained"]
    uniq_list = list(set(premature_stop_list))
    filtered_list = filter(None, uniq_list)
    print "Gene Models with premature stops\t%s" % len(filtered_list) 


    ##find the number of gene models with splice_donor or splice acceptor
    splice_don_accept = list(set(dataDict["splice_donor_variant"]) | set(dataDict["splice_acceptor_variant"]))
    uniq_list = list(set(splice_don_accept))
    filtered_list = filter(None, uniq_list)
    print "Gene Models with splice donor or acceptor variants\t%s" % len(filtered_list) 

    ### find the number of gene models in either splice_donor_variant, splice_acceptor_variant, or stop_gained
    stop_list = list(set(dataDict["splice_donor_variant"]) | set(dataDict["splice_acceptor_variant"]) | set(dataDict["stop_gained"]))
    uniq_list = list(set(stop_list))
    filtered_list = filter(None, uniq_list)
    uniqstop_list = filtered_list
    print "Gene Models in either splice donor/acceptor or stop_gained (combined truncations)\t%s" % len(filtered_list) 
    out_fh = open("truncation_list.tsv", "w")
    for model in sorted(filtered_list):
        out_fh.write(model + "\n")

    ### find the number of gene models in either splice_donor_variant, splice_acceptor_variant AND stop_gained
    stopandsplice_list = list((set(dataDict["splice_donor_variant"]) | set(dataDict["splice_acceptor_variant"])) & set(dataDict["stop_gained"]))
    uniq_list = list(set(stopandsplice_list))
    filtered_list = filter(None, uniq_list)
    print "Gene Models in either splice donor/acceptor AND stop_gained\t%s" % len(filtered_list)


    ## find the number of gene models in splice_region_variants
    splice_list = dataDict["splice_region_variant"]
    uniq_list = list(set(splice_list))
    filtered_list = filter(None, uniq_list)
    uniqsplice_list = filtered_list
    print "Gene Models with at least 1 splice_region_variant\t%s" % len(filtered_list) 

    ## additional gene models with splice_region_variants
    #addmodels_list = list(set(uniqstop_list) - set(uniqsplice_list))
    addmodels_list = list(set(uniqsplice_list) - set(uniqstop_list))
    uniq_list = list(set(addmodels_list))
    filtered_list = filter(None, uniq_list)
    print "Additional (compared to truncation list) Gene Models with 1 splice_region_variant\t%s" % len(filtered_list) 


    ## find the number of gene models with missense
    aa_list = dataDict["missense_variant"]
    uniq_list = list(set(aa_list))
    filtered_list = filter(None, uniq_list)
    uniqaa_list = filtered_list
    print "Gene Models with a missense\t%s" % len(filtered_list) 
    out_fh = open("missense_list.tsv", "w")
    for model in sorted(filtered_list):
        out_fh.write(model + "\n")


    truncs_not_missense = list(set(uniqstop_list) - set(uniqaa_list))
    print "Gene Models with truncations not overlapping with missense\t%s" % len(truncs_not_missense) 


    missense_not_truncs = list(set(uniqaa_list) - set(uniqstop_list) )
    print "Gene Models with missense and no truncations\t%s" % len(missense_not_truncs) 


    ### find the number of gene models with trunc and missense
    missenseandtrunc_list = list(set(uniqstop_list) & set(uniqaa_list))
    print "Gene Models with truncation and missense\t%s" % len(missenseandtrunc_list) 

    ##Find gene models not missense nor truncs
    notmissense_nor_truncs = list(set(stop_list) | set(aa_list))
    uniq_list = list(set(notmissense_nor_truncs))
    filtered_list = filter(None, uniq_list)
    notmissense_nor_truncs = list(set(gm_list) - set(filtered_list))
    print "Gene Models without missense nor truncation\t%s" % len(notmissense_nor_truncs) 

    
    ##include initiator codon
    init_list = dataDict["initiator_codon_variant"]
    uniq_list = list(set(init_list))
    filtered_list = filter(None, uniq_list)
    print "Gene Models with initiator_codon_variant\t%s" % len(filtered_list) 

    ##gene models that are downstream/upstream that don't have any other mutation
    downup_list = list()
    downup_list.extend(dataDict["downstream_gene_variant"])
    downup_list.extend(dataDict["upstream_gene_variant"])
    uniq_list = list(set(downup_list))
    filtered_list = filter(None, uniq_list)
    uniqdownup_list = filtered_list
    allothers_list = list()
    for conseq_type in dataDict:
        if conseq_type != "downstream_gene_variant" and conseq_type != "upstream_gene_variant":
            allothers_list.extend(dataDict[conseq_type])
    uniq_list = list(set(allothers_list))
    filtered_list = filter(None, uniq_list)
    uniqallothers_list = filtered_list
    final_list = list(set(uniqdownup_list) - set(uniqallothers_list))
    print "Gene Models with downstream/upstream but no other mutation tag\t%s" % len(final_list) 
