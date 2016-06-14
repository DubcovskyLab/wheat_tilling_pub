#!/usr/bin/env python

###############################################################
##Written by Tyson Howell 11-7-2015 UC Davis
#
#Version 1.6 updated 12-10-2015
#
#This is a script to rescue reads that have multiple mapping 
#locations which results in a low quality score. The output 
#matches the bwa sam format, but quality scores are adjusted
#to 255 (unknown quality) for rescued reads, and an additional
#tag XQ:i:# is added to indicate the original mapping quality
#
##NOTE:
##-This has only been tested on bams generated from bwa aln. It
##will not work with output from bwa MEM without modification.
##
##-This assumes that rescued mates (XT:A:M) are uniquely mapped
##
##-This first selects the contig with the smallest edit distance,
##followed by selecting contigs with the greatest number of
##matches/mismatches (based on the 'M' in CIGAR string) to avoid
##selecting mappings with an indel over a mapping with a mismatch
##
##-This deals with mappings one read at a time. Mate info
##may be incorrect after running this script!
##
##-This will not rescue mappings where there are 10 or more 
##"N" bases in the reference sequence. (XT:A:N flag)
##
##-When a read maps to multiple locations on a contig equally
##well, the read with the highest position is chosen. When 
##a read maps to two contigs of the same length, the contig
##that occurs last alphabetically is chosen
#
###############################################################

import argparse
import sys
import re

def set_bit(v, index, x):
	"""Set the index:th bit of v to x, and return the new value."""
	mask = 1 << index
	v &= ~mask
	if x:
		v |= mask
	return v

def index_containing_substring(the_list, substring):
	for i, s in enumerate(the_list):
		if substring in s:
			return i
	return -1

def revcomp(seq):
	complement = {
	'A': 'T',
	'B': 'V',
	'C': 'G',
	'D': 'H',
	'G': 'C',
	'H': 'D',
	'K': 'M',
	'M': 'K',
	'N': 'N',
	'R': 'Y',
	'S': 'S',
	'T': 'A',
	'V': 'B',
	'W': 'W',
	'Y': 'R',

	'a': 't',
	'b': 'v',
	'c': 'g',
	'd': 'h',
	'g': 'c',
	'h': 'd',
	'k': 'm',
	'm': 'k',
	'n': 'n',
	'r': 'y',
	's': 's',
	't': 'a',
	'v': 'b',
	'w': 'w',
	'y': 'r',
	}

	newseq = ''.join([complement[base] for base in list(seq)][::-1])
	return newseq

parser = argparse.ArgumentParser()
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-i", "--inputsam", required=True, help="Sam file to be processed (use \"<(samtools view input.bam)\" for bam files")
requiredNamed.add_argument("-l", "--lengths", required=True, help="Input file with lengths of fasta entries in form \"Length WHITESPACE Name\"")
parser.add_argument("-U", "--unique", action="store_true", help="Rescue only uniqely mapped reads (XT:A:U and XT:A:M flags) [F]")
parser.add_argument("-R", "--repetitive", action="store_true", help="Rescue only reads mapped to multiple best locations (XT:A:R flag) [F]")
parser.add_argument("-A", "--all", action="store_true", help="Rescue all mappings that have a low mapping quality and an XA flag (XT:A:U, XT:A:M and XT:A:R flags. Does not rescue mappings with XT:A:N flag [T]")
parser.add_argument("-I", "--IgnoreUniqueInfo", action="store_true", help="Ignores XT flags, treats all reads as if they were XT:A:R. Implies -A option [F]")
args = parser.parse_args()

if not (args.unique or args.repetitive or args.all):
	args.all = True
	sys.stderr.write("Output type not set, defaulting to rescue all reads\n")
	
if (args.unique or args.repetitive) and (args.IgnoreUniqueInfo):
	print("-U and -R are incompatible with the -I flag. Exiting...")
	sys.exit()
	
if (sum(map(bool, [args.unique, args.repetitive, args.all])) != 1):
	print("Please specify at most one of -U/--unique -R/--repetitive or -A/--all")
	sys.exit()
	
insam = open (args.inputsam)
lengthsfile = open (args.lengths)

length = {}
for line in lengthsfile:
	line = line.rstrip("\r\n")
	if line is '': continue #skip empty lines at end of file
	contig = line.split()
	length[contig[1]] = int(contig[0])

for line in insam:
	line = line.rstrip("\r\n")
	fields = line.split()
	#check if we want this line
	if ("XA:Z:" in line and int(fields[4]) < 20): #Select the reads we want to rescue
		#print "[READ]:", fields[0]
		#print line
		flags = {} #dictionary for storing relevant BWA flags
		mappings = [] #list for storing all identified mapping locations
		
		#parse additional fields
		for field in fields[11:]:
			flags[field[:2]] = field[5:]
		
		#check if the current mapping is optimal
		if (flags['XT'] is "U" or flags['XT'] is "M") and (args.all or args.unique) and not (args.IgnoreUniqueInfo): #Mapping is unique, just has quality <20 for some reason
			#Update quality to 255 and output
			fields.append("XQ:i:%s" % (fields[4])) #Add tag for real mapping quality
			fields[4] = str(255) #Change mapQ to 255
			#print fields[2] + " maps uniquely, looks good!"
			#print line
			print '\t'.join(fields)
			
		elif ((flags['XT'] is "R") and (args.all or args.repetitive)) or (args.IgnoreUniqueInfo and (flags['XT'] is "R" or flags['XT'] is "U" or flags['XT'] is "M")): #Read is mapped to multiple equally good locations
			
			#If there are more than 10 best locations, ignore it (this shouldn't happen since we only have the XA flag for up to 10 by default)
			if ('X0' in flags) and (int(flags['X0']) > 10): #There is no X0 field for XT:A:M alignments, so the extra test is needed to avoid a key error
				continue
			
			##Now we need to look at the alternative mapping locations
			#Put primary mapping in same format as alternates, build list of all mappings
			mainposition = None
			if (int(fields[1]) & 0x0010): #Read is mapped on reverse strand
				#print "Negative %s" % fields[1]
				mainposition = "-" + fields[3]
			else: #Read is mapped on forward strand
				#print "positive %s" % fields[1] 
				mainposition = "+" + fields[3]
			mappings = [[fields[2], mainposition, fields[5], flags["NM"]]]
			#0:chr, 1:pos, 2:CIGAR, 3:NM
			for info in flags['XA'].split(";"):
				if info is '': continue #There is a ';' after the last alternate hit which creates an empty element
				mappings.append(info.split(",")) #add alternate mappings

			#subset alternates by those that have the fewest mismatches
			minMM = min(int(MM[3]) for MM in mappings) #lowest number of mismatches (edit distance)
			minlist = [] #List of candidate alternative locations
			for el in mappings:
				if (int(el[3]) == minMM):
					minlist.append(el)
			#Need to build a new list to deal with corner cases of multiple locations on one contig, plus select the 'best' CIGAR strings 
			precandidates = []
			for contig in minlist:
				#Fields in contig 0:chr, 1:pos, 2:CIGAR, 3:NM
				CIGARmatches = re.findall(r'([0-9]+)([MIDNSHPX=])', contig[2]) #Parse cigar string into list of tuples [(number,matchtype), (number,matchtype)]
				matches = 0
				for group in CIGARmatches:
					if (group[1] is 'M' or group[1] is '='): matches += int(group[0])
				precandidates.append([length.get(contig[0]), contig[0], abs(int(contig[1])), contig[1], contig[2], contig [3], matches])
				#Fields in [pre]candidates 0:length(int) 1:contig name(string) 2: position(abs,int) 3:position(string) 4:CIGAR(string) 5:NM(string) 6:CIGAR number mismatches (int)
			
			#Filter candidates by number of matches
			maxmatches = max(contig[6] for contig in precandidates)
			candidates = []
			for el in precandidates:
				if (el[6] == maxmatches):
					candidates.append(el)
			candidates.sort(reverse=True) #Sorts by contig length, then by contig name (reverse alphabetically), then by position (reverse)
				
			#Check if the selected contig is the primary mapping
			if (fields[2] is candidates[0][1] and candidates[0][3] is mainposition):
				#Selected mapping is same as primary, update quality score and output
				fields.append("XQ:i:%s" % (fields[4])) #Add tag for real mapping quality
				fields[4] = str(255) #Change mapQ to 255
				
				#print fields[2] + " primary is already the best place to map, looks fine!"
				#print line
				print '\t'.join(fields)

				# Col   Field   Description
				# 0	 QNAME   Query (pair) NAME
				# 1	 FLAG	bitwise FLAG
				# 2	 RNAME   Reference sequence NAME
				# 3	 POS	 1-based leftmost POSition/coordinate of clipped sequence
				# 4	 MAPQ	MAPping Quality (Phred-scaled)
				# 5	 CIAGR   extended CIGAR string
				# 6	 MRNM	Mate Reference sequence NaMe ("=" if same as RNAME)
				# 7	 MPOS	1-based Mate POSistion
				# 8	 ISIZE   Inferred insert SIZE
				# 9	 SEQ	 query SEQuence on the same strand as the reference
				# 10	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
				# 11	OPT	 variable OPTional fields in the format TAG:VTYPE:VALUE

			else: #We need to update the primary mapping info
				#Update contig name
				if (fields[2] != candidates[0][1]): #If the contigs are different
					if (fields[6] is "="): #If the mate mapped to the same contig as the primary
						fields[6] = fields[2] #Change from "=" to old primary contig name
					fields[2] = candidates[0][1] #Change primary contig to the selected contig
					if (fields[2] is fields[6]): #If the new contig matches the alternate 
						fields[6] = "=" #update to "="
				
				#Update position
				fields[3] = str(candidates[0][2])
				
				#Update flag for strandedness
				if (candidates[0][3][:1] is "+"):
					#print "Positive: " + fields[1]
					fields[1] = str(set_bit(int(fields[1]), 4, 0))
					#print fields[1]
					
				else: #Strand is -
					#print "Negative: " + fields[1]
					fields[1] = str(set_bit(int(fields[1]), 4, 1))
					#print fields[1]
					
				#Update seq if needed
				if (mainposition[:1] != candidates[0][3][:1]): #If primary mapping is not the same strand as the candidate
					origseq = str(fields[9])
					fields[9] = revcomp(origseq) #reverse complement sequence
					origqual = str(fields[10]) #also reverse qualities
					fields[10] = origqual[::-1]
					#print origseq
					#print "New Seq: " + str(fields[9]) + " " + fields[10]
				#else:
					#print "Same Seq: " + str(fields[9]) + " " + fields[10]
					
				#Update mapQ
				fields.append("XQ:i:%s" % (fields[4])) #Add tag for real mapping quality
				fields[4] = str(255) #Change mapQ to 255
				#Update CIGAR
				fields[5] = candidates[0][4]
				#Update NM flag
				flags["NM"] = candidates[0][5]
				#Update XA flag
					#build list of new XA flags
				newXA = []
				for contig in mappings:
					#Fields in candidate 0:length(int) 1:contig name(string) 2: position(abs,int) 3:position(string) 4:CIGAR(string) 5:NM(string) 6:CIGAR number mismatches (int)
					#Fields in mappings 0:chr, 1:pos, 2:CIGAR, 3:NM
					#print contig
					#print [candidates[0][1], candidates[0][3], candidates[0][4], candidates[0][5]]
					if (contig != [candidates[0][1], candidates[0][3], candidates[0][4], candidates[0][5]]): #If the current contig is not the top candidate
						newXA.append(",".join(contig))
				XAindex = index_containing_substring(fields,"XA:Z:") #Get index of XA field (it is not always the same)
				fields[XAindex] = "XA:Z:" + ";".join(newXA) + ";" #Replace XA field. Append extra semicolon to match bwa sam format
				#print "Fixed this one!"
				#print line
				print '\t'.join(fields)
		elif(flags['XT'] is "N"): #These have 10 or more 'N's in the reference
			continue
		else:
			sys.stderr.write("No XT flag or incorrect XT flag, check input!\nOffending line:\n")
			sys.stderr.write(line)
			sys.stderr.write("\n")
			
