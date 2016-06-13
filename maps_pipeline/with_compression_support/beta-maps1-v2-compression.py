#! /usr/bin/env python2.6

import os, sys, math, gc, time
import threading, multiprocessing
from multiprocessing import Process, Queue
from optparse import OptionParser
from collections import defaultdict
import subprocess
import gzip

#Modified by TH 11-24-2015 for compressed input

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2012
# This work is the property of UC Davis Genome Center - Comai Lab
#
# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------
#
#MAPS - Part 1, maps-part1.py
#
#This program takes a parsed mpileup file and generates a "interesting" position file that serves as input to the second part of the MAPS pipeline. This first step is typically run once while the second program can be run quickly many times to try different mutation cutoffs.
#
#INPUT: This program takes a parsed mpileup.txt file as input.
#
#OUTPUT: 
#i. First, the user specifies an output file name for the positions of interest, this is used as input for the second MAPS program. 
#ii. There is an assay-[input file name] file generated that is an assay of the coverage of all "good" positions that pass the basic cutoffs. 
#iii. There is a snps-[input file name] file generated that records any "good" position that is a SNP position (a position that is polymporphic between the reference genome and the samples) across all libraries that have data for that position. 
#iv. There is a types-[input file name] file generated that outputs the counts of the different types of positions found. Only "good" positions are counted
#
#
#A position is defined as good if it passes the cutoffs for minimum position coverage, maximum position coverage, 
#and minimum of libraries with data for that position. This is also what we mean by assayed.
#
#What are the different types of positions: 
#i. het- These are "good" positions that pass the het cutoff parameters specified. (explained below) 
#ii. ref - These are "good" positions that match the reference. 
#iii. snp - These are positions that across all lib data for a position are a polymorphic between the reference and the samples. 
#iv. other - These are the "interesting" positions used in MAPS part 2. 
#v. rejected - These are positions that do not pass the basic position cutoffs mentioned to be a "good" position above. 
#vi. unknown - These are assayed positions for which the data that are too complex / strange to be useful in step 2.
#
#NOTE: 
#i. This first program reads the entire file into memory before running, so it is recommended not to run it on systems with limited memory. It typically requires > 1.5 times the size of the parsed mpileup file in RAM to run. For example, if the mpileup.txt file is 20 gbps, 30 Gbps of RAM are necessary to run the program on the whole file at once. Depending on the number of samples and the coverage, mpileup files can be very large and the required RAM will not be available. In this case, it is recommended to break the parsed mpileup file into smaller chunks to be processed separately, or leave them in chunks if they were already separated during the mpileup parsing step (i. e. by chromosome or scaffold. 
#ii. This program is threaded, so it can be used with the -t flag to specify the number of cores to be used.
#
#PARAMETERS, Default value in []: 
#1. REQUIRED: 
#i. -f or --file, The input parsed mpileup file. [required] 
#ii. -o or --out, The output interesting position destination file. [required] 
#2. OPTIONAL: 
#i. -t or --thread, Number of cores to be used while processing [1] 
#ii. -l or --MinLibs, Minimum number of libraries covered at least once to be considered a valid position. [3] 
#iii. --minCov or -v, Minimum total coverage to be considered as a valid position. Total coverage is the sum of the coverage at that position in all invididual libraries. [10] 
#iv. --maxCov or -c, Maximum total position coverage to be considered as a valid position. Total coverage is the sum of the coverage at that position in all invididual libraries. [2000] 
#v. --hetBothMinPer or -b, Minimum percentage for consideration of heterozygous position for both bases combined [95.00], i.e. % SNP base 1 and %SNP base 2 percent must sum to above this cutoff. This is based on the total coverage at a position. For example, if a base has the following coverage: 20A (50%), 16T (40%) and 4G (10%). SNP1 is A and SNP2 is T. SNP1 + SNP2 = 90%. 
#vi. --hetOneMinPer or -i, :Minimum percentage for consideration of heterozygous position for each of the two most common bases.[20.00] Each of the two SNP abse percentages must be at least this cutoff percent. 
#vii. --minCovNonMutBase or -u, Minimum coverage of a non-mutant base [5] 
#viii. --mode or -m, Output Mode: m == Mutation Detection, g == Genotyping [m]



usage = "\n\n%prog -f mpileup.txt -o outputfile.txt [OPTIONS] [-h for help]"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", dest="f", help="Input mpileup file.")
parser.add_option("-g", "--gff", dest="g",default = False, help="Input gff target file.")
parser.add_option("-o", "--out", dest="o", help="Output mpileup destination.")
parser.add_option("-t", "--threads", dest="numThreads", type = "int", default=1, help="Number of cores to use. (1)")
parser.add_option("-l", "--MinLibs", dest="minlibs", type = "int", default=3, help="Minimum number of libraries covered to be considered a valid position. (3)")
parser.add_option("--minCov", "-c", dest="mincov", type = "int", default=10, help="Minimum position coverage to be considered as a vaild position. (10)")
parser.add_option("--maxCov", "-C", dest="maxcov", type = "int", default=2000, help="Maximum position coverage to be considered as a vaild position. (2000)")
parser.add_option("--minSNPCut", "-n", dest="snpcut", type = "float",default=0.0, help="FOR TESTING USE: Minimum percent of total coverage for a snp to be considered relevant data. (0.0)")
parser.add_option("--hetBothMinPer", "-b", dest="hetBothMinPer", type = "float", default=95.00, help="Minimum percentage for consideration of het for both bases combined. (95.0)")
parser.add_option("--hetOneMinPer", "-i", dest="hetOneMinPer", type = "float", default=20.00, help="Minimum percentage for consideration of het for each of the two bases. (20.0)")
parser.add_option("--minCovNonMutBase", "-u", dest="minNonMut", type = "int", default=5, help="Minimum coverage of a nonmutant base. (5)")
parser.add_option("--mode", "-m", dest="mode",  type = "str", default='m', help="Output Mode:  m==Mutation Detection,  g==Genotyping (m) ")
parser.add_option("--mutHet", "-H", dest="hetout", action="store_true", default = False, help="In mutation mode output hets (Off)")
parser.add_option("-L", "--length", dest="in_len", help="File containing length of parsed mpileup in lines")

(opt, args) = parser.parse_args()
#open parsed mpileup input file

#bring in all command line options
numThreads = opt.numThreads
minlibs = opt.minlibs
maxcov = opt.maxcov
mincov = opt.mincov
snpcut = opt.snpcut
hetBothMinPer = opt.hetBothMinPer
hetOneMinPer = opt.hetOneMinPer
minNonMut = opt.minNonMut
mode = str(opt.mode)
args = [minlibs, maxcov, mincov, snpcut, hetBothMinPer, hetOneMinPer, minNonMut]

###############
#pieces to break oither into
chunks = 3 * int(round(256/opt.numThreads,0))
#chunks = 6


# Uses wc to get te number of lines in the file
def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

#split to length
def splitter(l, n):
    i = 0
    chunk = l[:n]
    while chunk:
        yield chunk
        i += n
        chunk = l[i:i+n]

#text formatting
def form(flo):
   return str(flo).split('.')[0]+'.'+str(flo).split('.')[1][:2]

#-------------------------------

#each thread takes 1/#threads porton of the file, using a queue
#to return the relevant data to preserve order without using temporary files
class MyThread (multiprocessing.Process):

   def __init__ (self, filen, res, startpos, endpos, args, vlist, chunks):
      self.q = Queue()
      self.r = Queue()
      self.filen = filen
      self.res = res
      self.startpos = startpos
      self.endpos = endpos
      self.args = args
      self.vlist = vlist
      self.chunks = chunks
      multiprocessing.Process.__init__ (self)     

   def run (self):
      #take in parameters
      minlibs, maxcov, mincov, snpcut, hetBothMinPer, hetOneMinPer, minNonMut = self.args  
      track = 0
      #all = len(self.lines)
      assay = {}
      targetassay = {}
      #setup empty assay position table, 0-maxcov per cov point
      for covnum in range(0,maxcov+1):
         assay[covnum] = {}
         targetassay[covnum] = {}
         for libname in liblist+['All']:
            targetassay[covnum][libname] = 0
            assay[covnum][libname] = 0
    
      counters = {'het': 0, 'ref': 0, 'snp': 0, 'other': 0, 'rejected': 0, 'unknown': 0}
      snps = []
      others = []
      het = []
      #for each position
      counter = 1
      ctgood = 0
      print self.startpos, self.endpos
      all = self.endpos - self.startpos
      f = gzip.open(self.filen)
      f.readline()
      for l in f:
         if counter < self.startpos or counter > self.endpos:
            counter +=1
            continue
            
#         if ctgood % 100008 == 0:
#            print self.res, ctgood,'/',all
#         ctgood +=1
         counter +=1

#      for l in self.lines:
         if track % 250008 == 8:
            print str(self.res)+'-\t'+str(track)+'/'+str(all), counters, sum(counters.values()) == track
         track +=1         
         xline = l[:-1].split('\t')
         chrom = xline[0]
         pos = xline[1]
         ref = xline[2]
         #splitting by lib
         check = list(splitter(xline[3:],4))
         # to hold 'valid' data for each lib
         newVals = []
         totcov = 0
         coveredlibs = 0
         firstSNP = []
         refCheck = 0
         libindex = 0
         covs = {}
         #set empty dict to hold cov by lib
         for x in liblist:
            covs[x] = 0
         #for the data in each lib (snp1, snp2, snp3, cov)
         for lib in check:
            templib = []
            #if empty data
            if lib[0] == '.':
               newVals.append(lib)
               templib = ['.','.','.','.']
            #else if not empty
            else:
               #FOR SNP FIXING<--- NOT REALLY USED -- -- - - -- -
               ocov = int(lib[-1])
               tempcov = ocov+0
               for snp in reversed(lib[:-1]):
                  if snp == '.':
                     templib = ['.'] + templib
                  else:
                     snpd = snp.split('_')
                     base = snpd[0]
                     bper = snpd[1]
                     
                     scov = int(round((float(bper)/100.00)*ocov))
                     sper = (scov/float(tempcov))*100.00
                    
                     if sper <= snpcut:
                        templib = ['.'] + templib
                        tempcov -= scov
                     else:
                        templib = [base+'_'+form(sper)] + templib
               #end snpcut - - -- --- - - ------------------------
               templib.append(str(tempcov))                  
               newVals.append(templib)
               totcov += int(templib[-1])
               coveredlibs+=1
               covs[liblist[libindex]] = templib[-1]
               firstSNP.append(templib[0])
               #if not ref, break all ref condition check
               if templib[0] !=  ref+'_100.0':
                  refCheck = 1
            libindex+=1
         firstSNP = list(set(firstSNP))
         #To be assayed needs >= mincov per pos, <= maxcov, and x # minlibs or >              
         #print "DEBUG:", totcov, coveredlibs
         if totcov <= maxcov and coveredlibs >= minlibs and totcov >= mincov:
            #for each lib, assay the cov at that position
            if opt.g == False:
               for lname in liblist:
                  val = int(covs[lname])
                  assay[val][lname] += 1
               assay[totcov]['All'] +=1
            else:
               key1 = int(pos) / 10000
               key2 = int(pos) % 10000
               try:
                  if key1 in self.vlist[chrom]:
                     check = key2 in self.vlist[chrom][key1]
                  else:
                     check = False
               except KeyError:
                  check = False
               #print check
               #print covs
               for lname in liblist:
                  val = int(covs[lname])
                  if check == False:                     
                     assay[val][lname] += 1                   
                  else:
                     targetassay[val][lname] += 1
                  if check == False:
                     assay[totcov]['All'] +=1
                  else:
                     targetassay[totcov]['All'] +=1                   
            #IF ALL THE SAME ----------------
            #IF REF +++++++++++
            if refCheck == 0:
               counters['ref'] += 1
            #IF SNP +++++++++++++
            elif len(firstSNP) == 1:
               counters['snp'] += 1
               line = '\t'.join([chrom,pos,ref]+firstSNP) +'\n'  
               snps.append(line)
            #IF DIFFERENT ----------------
            else:
               summary = {'A': [0, 0], 'C': [0, 0], 'G': [0, 0], '+': [0, 0], '*': [0, 0], 'T': [0, 0]}
               for mlib in newVals:
                  if mlib != ['.','.','.','.']:
                     baseCov = int(mlib[-1])
                     for snp in mlib[:-1]:
                        if snp != '.':
                           snpx = snp.split('_')
                           snpbase = snpx[0]
                           if '+' in snpbase or '-' in snpbase:
                              snpbase = '+'
                           snpPer = float(snpx[-1])/100.00
                           summary[snpbase][0]+=1
                           summary[snpbase][1] += int(round(float(snpPer)*baseCov))
               #now after summary generated
               #print summary      
               #check for het
               rsumm = [[z,summary[z][1]] for z in summary]
               rsumm.sort(lambda x1, y1: cmp(int(y1[1]), int(x1[1])))
               hetcheck = rsumm[0][1]+  rsumm[1][1]
               hetkeep = rsumm
               top2 = [x[0] for x in rsumm[:2]]
               libcounts = [summary[x][0] for x in top2]
               hetter = 0
               #both have to be in more than 1?
               if len(filter(lambda x: x>1, libcounts)) == 2:
                  hetter = 1
               #IF HET +++++++++++              
               if hetcheck/float(totcov)*100 >= hetBothMinPer and rsumm[0][1]/float(totcov)*100 >= hetOneMinPer and rsumm[1][1]/float(totcov)*100 >= hetOneMinPer and hetter == 1:
                     if mode == 'm':
                        counters['het'] +=1
                        if opt.hetout == True:
                           line = [chrom, pos, ref]
                           line += [totcov]+ summary['A'] + summary['T'] + summary['C'] + summary['G'] + summary['*']+ summary['+']
                           for sub in newVals:
                              if sub[1] == '.':
                                 line+= [sub[0],sub[-1]]
                              else:
                                 line+=[sub[0]+'-'+sub[1],sub[-1]]
                           line = map(lambda x: str(x), line)
                           het.append('\t'.join(line)+'\n')
                     elif mode == 'g':
                        counters['other'] +=1
                        line = [chrom, pos, ref]
                        line += [totcov]+ summary['A'] + summary['T'] + summary['C'] + summary['G'] + summary['*']+ summary['+']
                        for sub in newVals:
                           if sub[1] == '.':
                              line+= [sub[0],sub[-1]]
                           else:
                              line+=[sub[0]+'-'+sub[1],sub[-1]]
                        line = map(lambda x: str(x), line)
                        others.append('\t'.join(line)+'\n')
                     else:
                        print "What????"
                        parser.error("Het: Please check your command line paramters with -h or --help")                      
               #EVERYTHING ELSE ++++++++++ 
               else:      
                  unkno = 0
                  if 1 in libcounts:
#                     onelibcov = filter(lambda x: summary[x[0]][0]==1, rsumm)[0][1]
#                     combinedcov = sum([summary[x][1] for x in top2])
#                     checkval = combinedcov - onelibcov
#                     if checkval != rsumm[0][1]:
#                        print xline, '\n', onelibcov, '\n', combinedcov, '\n', checkval,'\n', summary,'\n', top2, '\n', rsumm, '\n', libcounts, '\n---------------------------'
                     if rsumm[0][1] < minNonMut:
                        counters['unknown'] +=1
                        unkno = 1                  
                  if unkno == 0:     
                     counters['other'] +=1
                     line = [chrom, pos, ref]
                     line += [totcov]+ summary['A'] + summary['T'] + summary['C'] + summary['G'] + summary['*']+ summary['+']
                     for sub in newVals:
                        if sub[1] == '.':
                           line+= [sub[0],sub[-1]]
                        else:
                           line+=[sub[0]+'-'+sub[1],sub[-1]]
                     line = map(lambda x: str(x), line)
                     others.append('\t'.join(line)+'\n')  
         else:
            counters['rejected'] += 1
      print str(self.res)+" Put Other"
      
      tempct = len(others)/self.chunks
      for iter1 in range(self.chunks):
         print str(self.res)+" Put Other "+str(iter1+1)
         if iter1 == self.chunks-1:
            self.q.put(others[iter1*tempct:])
         else:
            self.q.put(others[iter1*tempct:(iter1+1)*tempct])
         checkit = self.r.get()

      
      
      
      print str(self.res)+" Put Rest", checkit

      self.q.put(het)
      self.q.put(snps)
      self.q.put(assay)
      self.q.put(targetassay)
      self.q.put(counters)
      print str(self.res)+" Done"

#def getkeyset(pos):
#   #5000
#   key1 = pos / 10000
#   rem = pos % 10000
#   key2 = 


try:
   f = gzip.open(opt.f)
   filename = opt.f
except:    
   parser.error("f: Please check your command line paramters with -h or --help") 

#open output file 
validlist = {}
try:
   o = open(opt.o, 'w')
   nameo = opt.o  
except:    
   parser.error("o: Please check your command line paramters with -h or --help")

bct = 0
chromlist = []
if opt.g != False:
   print "Reading target GFF"
   g = open(opt.g)
   for l in g:
      x = l.split('\t')
      chrom = x[0]
      if chrom not in chromlist:
         print chrom
         chromlist.append(chrom)
         validlist[chrom] = {}
      spos = x[3]
      epos = x[4]
      valid = range(int(spos), int(epos)+1)
      for act in valid:
         bct +=1
         key1 = act / 10000
         key2 = act % 10000
         if key1 not in validlist[chrom]:
            validlist[chrom][key1] = []
         validlist[chrom][key1].append(key2)

      


   
#Make headers ------------------
head = f.readline()
totlibs = len(head.split('\t')[3:])/4
head = head[:-1].split('\t')
newhead = ["Chrom/Scaffold", "Pos", "Ref"]
#liblist = map(lambda x: x.split('-')[-1],head[3::4])
liblist = map(lambda x: ''.join(x.split('-')[1:]),head[3::4])
newhead += ['TotCov','A-#Lib','A-Cov', 'T-#Lib','T-Cov', 'C-#Lib','C-Cov', 'G-#Lib','G-Cov', '*-#Lib','*-Cov', '+-#Lib','+-Cov']

temp = list(splitter(head[3:],4))
for x in temp:
   newhead += [x[0],x[-1]]
#-----------

snphead = '\t'.join(["Chrom/Scaffold", "Pos", "Ref", "SNP"])+'\n'
outhead = '\t'.join(newhead)+'\n'


f.close()

flen = None
if opt.in_len:
    with open(opt.in_len, 'r') as l:
        print l
        flen = int(l.readline().strip().split()[0])
else:
    flen = file_len(filename)
#flen = 340046839
cutnum = (flen-1)/numThreads+1
cutset = []
for i in range(numThreads):
   cutset.append([i*cutnum+1, min((i+1)*cutnum, flen)])



#read in file
#a = []
#for l in f:
#   a.append(l)
#split file in to #ofCores segments, preserving order

#a = list(splitter(a,len(a)/numThreads+1))
##del(a)
#gc.collect()

countert = 0
threads = {}
print cutset
for x in cutset:
   countert+=1
   print countert
   threads[countert] = MyThread(filename, countert, x[0], x[1], args, validlist, chunks)
   threads[countert].start() 
#spin out the threads
#for x in a:
#   counter+=1
#   print counter
#   threads[counter] = MyThread(x, counter, args)
#   threads[counter].start() 
   
#f.close()
#del(a)
#gc.collect()

#set lists for results
all = []
allothers = []
allsnps = []
allhets = []
allassay = []
alltargetassay = []
allcounters = []
#using a queue, get all values then join thread,
#once again in order
for x in range(1, countert+1):
   print "Getting others ", x
   for inter2 in range(chunks):
      allothers+=threads[x].q.get()
      threads[x].r.put("continue")
   print "Getting rest ", x
   allhets+=threads[x].q.get()
   allsnps+=threads[x].q.get()
   allassay.append(threads[x].q.get())
   alltargetassay.append(threads[x].q.get())
   allcounters.append(threads[x].q.get())
   threads[x].join()
   print x, "joined"
 
#write others lines -------------------------------------------------------
o.write(outhead)
for x in allothers:
   o.write(x)
o.close()

# write snps --------------------------------------------------------------

if opt.f.endswith('.gz'):
    checkedname = os.path.splitext(opt.f)[0]
else:
    checkedname = opt.f

try:
   t = open("snps-"+checkedname,'w')
except:    
   parser.error("SNPs: Please check your command line paramters with -h or --help")

t.write(snphead)
for x in allsnps:
   t.write(x)

t.close()

if opt.hetout == True:
   # write hets --------------------------------------------------------------
   try:
      ts = open("hets-"+checkedname,'w')
   except:    
      parser.error("hetss: Please check your command line paramters with -h or --help")
   
   ts.write(outhead)
   for x in allhets:
      ts.write(x)
   
   ts.close()

# write out type counts -----------------------------------------------------
try:
   d = open("type-"+checkedname,'w')
except:    
   parser.error("type: Please check your command line paramters with -h or --help")
   
sumcounters = {'het': 0, 'ref': 0, 'snp': 0, 'other': 0, 'rejected': 0, 'unknown': 0}
for subset3 in allcounters:
   for label in sumcounters.keys():
      sumcounters[label] += subset3[label]

print sumcounters      
d.write('\n'.join([x+':\t'+str(sumcounters[x]) for x in sumcounters])+'\n')
d.close()


#join all assays and write out master assay ----------------------------------
try:
   if opt.g == False:
      j = open("assay-"+checkedname,'w')
   else:
      j = open("non-target-assay-"+checkedname,'w')
except:    
   parser.error("assay: Please check your command line paramters with -h or --help")

sumassay = {}
for covnum in range(0,maxcov+1):
   sumassay[covnum] = {}
   for libname in liblist+['All']:
      sumassay[covnum][libname] = 0

for set in allassay:   
   for cov in range(0,maxcov+1):
      for libname in liblist+['All']:
         sumassay[cov][libname] += set[cov][libname]
   
j.write('\t'.join(["Cov"]+liblist+['All'])+'\n')    
for val in range(1,maxcov+1):
   line = [str(val)]
   for ln in liblist+['All']:
      line.append(str(sumassay[val][ln]))
   j.write('\t'.join(line)+'\n')
j.close() 

#join all assays and write out master assay ----------------------------------
if opt.g != False:
   try:
      j = open("target-assay-"+checkedname,'w')
   except:    
      parser.error("targetassay: Please check your command line paramters with -h or --help")
   
   tsumassay = {}
   for tcovnum in range(0,maxcov+1):
      tsumassay[tcovnum] = {}
      for libname in liblist+['All']:
         tsumassay[tcovnum][libname] = 0
   
   for set in alltargetassay:   
      for cov in range(0,maxcov+1):
         for libname in liblist+['All']:
            tsumassay[cov][libname] += set[cov][libname]
      
   j.write('\t'.join(["Cov"]+liblist+['All'])+'\n')    
   for val in range(1,maxcov+1):
      line = [str(val)]
      for ln in liblist+['All']:
         line.append(str(tsumassay[val][ln]))
      j.write('\t'.join(line)+'\n')
   j.close() 
