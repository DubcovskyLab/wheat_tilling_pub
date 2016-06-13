#! /usr/bin/env python2.6

import os, sys, math
from optparse import OptionParser
import threading, multiprocessing
from multiprocessing import Process, Queue
import gzip

#Edited by TH 11-24-2015. Output is to stdout so that the data can be compressed on the fly, or piped. Intermidate files are also compressed to reduce I/O

#Ucdavis Genome Center
#Meric Lieberman, 2014
# This work is the property of UC Davis Genome Center

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------

#Part 1: run-mpileup.py

#This program is meant to be run on a directory of sorted.bam files with the corresponding .bai files. It will generate a mpileup file with columns for each library.

#INPUT:
#This program is run in a folder full of .sorted.bam files as input

#OUTPUT:
#This program outputs a mpileup file.

#NOTE:
#If the program samtools is not in /usr/bin, then the path to samtools must be specified using the command line parameters

#PARAMETERS, default value in []:
#1. REQUIRED:
#-r or reference_file, The alignment reference (fasta format) [required]
#-o or--output_file, The output mpileup.txt filename [required]
#2. OPTIONAL:
#-q or --mapqual, Minimum mapping quality for an alignment to be used [20]
#-Q or --basequal, Minimum base quality for a base to be considered [20]
#-d or --maxdepth, Max per-BAM depth coverage to avoid excessive memory usage [8000]
#-s or --samtools, File path to Samtools [/usr/bin/samtools]
#-t or --thread, Number of cores to be used while processing [1]

usage = "\n\n%prog -r reference.fa -o output.txt [-q y] [-Q x] [-s path to Samtools} > out.txt"
usage += "\nRun in a directory only full of sorted.bam files, will generate a mplieup from all bams. Outputs to STDOUT"
parser = OptionParser(usage=usage)
parser.add_option("-r", "--reference_file", dest="ref", help="Alignment reference file.")
parser.add_option("-o", "--output_file", dest="dest", help="Output file name.")
parser.add_option("-q", "--mapqual", dest="mapqual", default="20", help="(OPTIONAL, default = 20) Minimum mapping quality for an alignment to be used")
parser.add_option("-Q", "--basequal", dest="basequal", default="20", help="(OPTIONAL, default = 20) Minimum base quality for a base to be considered")
parser.add_option("-d", "--maxdepth", dest="maxdepth", default="8000", help="(OPTIONAL, default = 8000) Max per-BAM depth to avoid excessive memory usage")
parser.add_option("--samtools", "-s", dest="pathSAM",  type = "str", default='/share/apps/samtools-github-1.18/samtools', help="File path to Samtools")
parser.add_option("-t", "--threads", dest="threads", type = "int", default=1, help="Number of cores to use. (1)")
parser.add_option("--bamname", "-n", dest="bamFileName",  type = "str", default='_aln.sorted.bam', help="This is the file designation that the bam files must end with, the .bai files must match this convention as well default = _aln.sorted.bam")

(opt, args) = parser.parse_args()
mapqual = opt.mapqual
basequal = opt.basequal
maxdepth = opt.maxdepth



class MyThread (multiprocessing.Process):

   def __init__ (self, rline, idn):
      #quey out
      self.q = Queue()
      #response back for handshake
      self.r = Queue()
      #line to run
      self.runline = rline
      #thread identifier number
      self.idn = idn
      multiprocessing.Process.__init__ (self)     

   def run (self):
      sys.stderr.write("Thread "+str(self.idn)+" running\n")
      sys.stderr.write("%s\n" % self.runline)
      os.system(self.runline)
      #handshake for done
      self.q.put("Done")
      check = self.r.get()



def readFastaLen(filename):
   f = open(filename)       
   #read in a fasta
   name = ''
   seq = 0
   seqs = {} 
   while True:
      l = f.readline()
      if l == '':
         seqs[name] = seq
         break
      elif l[0] == '>':
         if name != '':
            seqs[name] = seq
         seq = 0
         name = l.replace('\n','').replace('>','')
      else:
         seq+=len(l.replace('\n',''))
   f.close()
   return seqs

###############
#Partition ref seqs
#########

#Getting all ref sizes
try:
   sizes = readFastaLen(opt.ref)
except:
   parser.error("Could not read reference. Please check your command line paramters with -h or --help")

#
#for x in all:
#   sizes[x] = len(all[x])

#figuring out size cutoff
clist = sizes.keys()
clist.sort()
allvals = sum(sizes.values())
piecelen = allvals/opt.threads

#breaking into #of threads chinks
varlist = [[]]
temp = 0
for c in clist:
   tn = sizes[c]
   #For last block (may be a bit bigger to pick up slack)
   if len(varlist) == opt.threads:
      varlist[-1].append(c)
      temp+=tn      
   #If bigger than cutoff, close this block and start nw
   elif temp + tn < piecelen*1.1:
      varlist[-1].append(c)
      temp+=tn
   #else, still smaller, add to existing block
   else:
      varlist+= [[c]]
      temp = tn
#For size testing
#fsize = []
#for x in varlist:
#   fsize.append(sum(map(lambda z: sizes[z], x)))
   
#Create temp file name standards 
front = "temp-region-"
back = ".bed"
ofront = "temp-mpileup-part-"
oback = ".txt"
#regions are the temprary text .bed files
#outputs are the result fractional mpileup files
regions = []
outputs = []
for x in range(1, len(varlist)+1):
   treg = front+str(x)+back
   tdes = ofront+str(x)+oback
   regions.append(treg)
   outputs.append(tdes)
   ot = gzip.open(treg, 'w')
   #Write out refname\t1\tlen(ref)
   for item in varlist[x-1]:
      line = [item, '1', str(sizes[item])]
      ot.write('\t'.join(line)+'\n')
   ot.close()


#check to make sure input and output files entered correctly
file = opt.ref
if opt.dest != '-':
    try:
       file = opt.ref
       o = open(opt.dest, 'w')
    except:
       parser.error("Could not open output. Please check your command line paramters with -h or --help")

li = os.listdir(os.getcwd())
ind = filter(lambda x: x.endswith(opt.bamFileName), li)
ind.sort()
a = map(lambda x: ["Cov-"+x.split('_')[0].replace('lib','').replace('Lib',''),'Call-'+x.split('_')[0].replace('lib','').replace('Lib',''),'Qual-'+x.split('_')[0].replace('lib','').replace('Lib','')], ind)
b = [item for sublist in a for item in sublist]
header = '\t'.join(['Chrom', 'Pos', 'Ref']+b)+'\n'


#write out header for final combined file
if opt.dest == '-':
    sys.stdout.write("%s" % header)
else:
    o.write(header)
    o.close()

#generate individual run lines for each chunk
runlinebams = ' '.join(ind)
runlines = []
for x in zip(regions, outputs):
   #line = opt.pathSAM + " mpileup -l "+x[0]+" -d "+ maxdepth +" -Q "+basequal+" -q "+mapqual+" -f "+file+" "+runlinebams+" > "+x[1]
   line = opt.pathSAM + " mpileup -l "+x[0]+" -d "+ maxdepth +" -Q "+basequal+" -q "+mapqual+" -f "+file+" "+runlinebams+" | gzip > "+x[1]
   runlines.append(line)

#spin out threads from 1 to numthreads
countert = 0
threads = {}
for x in runlines:
   countert+=1
   sys.stderr.write("Thread %s\n" % countert)
   threads[countert] = MyThread(x,countert)
   threads[countert].start() 


#using a queue, get all values then join thread,
#once again in order
#uses a put() and get() in each direction to handshake
for x in range(1, countert+1):
   sys.stderr.write( "Checking Sub %s\n" % str(x))
   check = threads[x].q.get()
   threads[x].r.put("continue")
   #print "Ending Sub ", x
   threads[x].join()
   sys.stderr.write( "Sub %s joined.\n" % str(x))

#all children have been gathered. cat together IN ORDER
if opt.dest == '-':
    catline = "zcat "+' '.join(outputs)
    sys.stderr.write("%s\n" % catline)
    os.system(catline)
else:
    catline = "zcat "+' '.join(outputs)+ " >> "+opt.dest
    print catline
    os.system(catline)

#clean up intermediary files
os.system("rm temp-region-*.bed")
os.system("rm temp-mpileup-part-*.txt")

