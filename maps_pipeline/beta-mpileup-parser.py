#! /usr/bin/env python2.6

import os, sys, math, datetime, gc, time
import threading, multiprocessing
from optparse import OptionParser
import subprocess

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2012
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------

#Part 2: mpileup-parser.py
#
#This program parses a mpileup file to a simplified format to be used with our MAPS 
#mutation and genotyping package.
#
#INPUT:
#This program takes an mpileup.txt file as input
#
#OUTPUT:
#This program outputs a parsed mpileup.txt file.
#
#NOTE:
#This program reads the entire file into memory before parsing, so it is recommended not to be run on systems with limited memory. It typically requires 1.5 times the size of the mpileup file in RAM to run. Many machines will not have this, and in this case it is recommended to break the mpileup file into smaller chunks to be processed separately. (i. e. by chromosome or scaffold)
#If being used with the MAPS package: the first step in MAPS is also threaded so it may be best to leave the chunks separate and process them individually in MAPS as well. Results can be combined at the end without compromising the results.
#This program is threaded, so it can be used with the -t flag to specify the number of cores to be used
#
#PARAMETERS:
#1. REQUIRED, default value in []:
#-f or --mpileup_file, The input mpileup file. [required]
#2. OPTIONAL:
#-t or --thread, Number of cores to be used while processing. [1]

start = time.time()

usage = "\nUSAGE: %prog [-t #threads] -f mpileup_file.txt"
parser = OptionParser(usage=usage)
parser.add_option("-t", "--thread", dest="threads", default="1", help="How many threads to use during processing. (DEFAULT == 1")
parser.add_option("-f", "--mpileup_file", dest="file", help="Input mpileup file file.")
(opt, args) = parser.parse_args()

numThreads = int(opt.threads)

path = "/share/scripts/"
path = "./"

#split file into chunks
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
#accepted base
def test(s):
   valid = ['a','t','c','g','A','T','C','G','.',',','*','n','N']
   test = 0
   for x in valid:
      if x+'+' in s:
         test = 1
         break
   if test == 1:  
      return 1
   else:
      return 0

#parse one pileup line
#based on parsing script from Joeseph Fass

class MyThread (multiprocessing.Process):

   def __init__ (self, filen, res, startpos, endpos):
      self.filen = filen
      self.res = res
      self.startpos = startpos
      self.endpos = endpos
      multiprocessing.Process.__init__ (self)
      
   def run (self):
      counter = 1
      ctgood = 0
      print self.startpos, self.endpos
      all = self.endpos - self.startpos
      ot = open("temp-parse-"+str(self.res)+".txt",'w')
      f = open(self.filen)
      f.readline()
      for l in f:
         if counter < self.startpos or counter > self.endpos:
            counter +=1
            continue
            
         if ctgood % 100008 == 0:
            print self.res, ctgood,'/',all
         ctgood +=1
         counter +=1
         #doOneLine(self.lines[k], ot) 
         #y = self.lines[c]
         ######
         #def doOneLine(y, fh):
         k = l[:-1].split('\t')
         refseq = k[0]
         position = k[1]
         refbase = k[2]
         div = list(splitter(k,3))
         result = div[0]
         #for a lib
         for sub in div[1:]:
            depth = sub[0]
            changes = sub[1]
            qualities = sub[2]
            try:
               mean_SQ = form(sum(map(lambda x: ord(x)-33, list(qualities)))/float(len(qualities)))
            except ZeroDivisionError:
               mean_SQ = "0.0"
            if float(mean_SQ) < 20.00:
               result += ['.','.','.','.']
               continue
            else:  
               inserts = {}
               quals = {'a':0,'A':0,'c':0,'C':0,'g':0,'G':0,'t':0,'T':0,'.':0,',':0,'*':0}
               valid = ['a','A','c','C','g','G','t','T','.',',','*']
            
               #depth = t[3]
               #changes = t[4]
               #qualities = t[5][:]
               #mappingqual = t[6][:-1]
               count = 0
               index = 0
               x = 0
               temp = sub[1]
               temp = temp.replace('^+','')  
               while test(temp) == 1:
                  pin = temp.index('+')
                  numb = ""
                  i = 0
                  while 1:
                     if temp[pin+1+i].isdigit():
                        numb+= temp[pin+1+i]
                        i+=1
                     else:
                        break    
                  
                  #take = int(temp[pin+1])
                  take = int(numb)
                  #print take
                  total = temp[pin-1:pin+2+take]
                  cleantotal = '.'+str.upper(total[1:])
                  try:
                     inserts[cleantotal] +=1
                  except:
                     inserts[cleantotal] = 1
                     
                  temp = temp.replace(total, '')
                  sub[1] = sub[1].replace(total, '')
               if '+' in sub[1]:
                  pin = sub[1].index('+')
                  if sub[1][pin-1] != '^':
                     print "error, +/- found"
                     print sub[1]
            
               while 1: 
                  if x >= len(sub[1]):
                     break
                  elif sub[1][x] in valid:
                     quals[sub[1][x]] +=1
                     index+=1
                  elif sub[1][x] == "^":
                     x+=1
                  elif sub[1][x] == '+' or sub[1][x] == '-':  
                     temp = ""
                     i = 0
                     while 1:
                        if sub[1][x+1+i].isdigit():
                           temp+= sub[1][x+1+i]
                           i+=1
                        else:
                           break
                     x+= int(temp)+len(temp)
            
                  x+=1
            
               total_HQ = sum(quals.values())
               Aa_HQ = quals['A']+quals['a']
               Tt_HQ = quals['T']+quals['t']
               Cc_HQ = quals['C']+quals['c']
               Gg_HQ = quals['G']+quals['g']
               match_HQ = quals['.']+quals[',']
               dels = quals['*']
               
         
               try:
                  if len(inserts) >0:
                     #print mean_SQ, qualities
                     inlist =  map(lambda x: [x, inserts[x]], inserts.keys())
                     inlist.sort(lambda x, y: cmp(y[1], x[1]))
                     inname = inlist[0][0]
                     incount = inlist[0][1]  
                     inper = incount/float(incount + total_HQ)*100
                     delper = dels/float(incount + total_HQ)*100
                     oin = str(sum(map(lambda x: x[1], inlist)))
                     total_HQ += int(oin)
                     scan = [[inname,inper]]
                  else:
                     inname = '.' 
                     inper = '.'
                     oin = '.'
                     delper = dels/float(total_HQ)*100
                     scan = []
              
                  aPer = Aa_HQ/float(total_HQ)*100
                  tPer = Tt_HQ/float(total_HQ)*100
                  cPer = Cc_HQ/float(total_HQ)*100
                  gPer = Gg_HQ/float(total_HQ)*100
                  matchPer = match_HQ/float(total_HQ)*100
                  scan += [['A',aPer],['T',tPer],['C',cPer],['G',gPer],['*',delper]]
                  for w in scan:
                     if w[0] == refbase:
                        w[1] += matchPer
                  scan.sort(lambda x, y: cmp(float(y[1]),float(x[1])))
                  scan2 = []
                  for w in scan:
                     if w[1] == 0.0:
                        scan2.append('.')
                     else:
                        scan2.append(w[0]+'_'+form(w[1]))
                  
                  result+= [scan2[0], scan2[1], scan2[2], str(total_HQ)]
               except ZeroDivisionError:
                  result+=['.','.','.','.']
      
         ot.write('\t'.join(result)+'\n')         
         #######    
      ot.close()
      f.close()

      
# Uses wc to get te number of lines in the file
def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])
    
    
     



      
t1 =  datetime.datetime.now()      
#print t1 
try: 
   filename = opt.file  
   f = open(filename)
   o = open("parsed_"+filename.split('/')[-1],'w')
except:
   parser.error("Please check your command line paramters with -h or --help")
flen = file_len(filename)
cutnum = (flen-1)/numThreads+1
cutset = []
for i in range(numThreads):
   cutset.append([i*cutnum+1, min((i+1)*cutnum, flen)])
 
t = f.readline()
if t == '':
   sys.exit("Empty pileup file")
f.seek(0)
header = f.readline()
header =header[:-1].split('\t')
h2 = list(splitter(header[3:],3))
newhead = header[:3]
libs = []
for lab in h2:
   lname = ('-'.join(lab[0].split('-')[1:]))
   libs.append(lname)
   newhead += ["SNP1-"+lname,"SNP2-"+lname,"SNP3-"+lname,"Cov-"+lname]

o.write('\t'.join(newhead)+'\n')
o.close()
f.close()
#
#a = []
#for l in f:
#   a.append(l)
#
#
#gc.collect()
counter = 0
threads = {}
cat = "cat"
print cutset
for x in cutset:
   counter+=1
   print counter
   threads[counter] = MyThread(filename, counter, x[0], x[1])
   cat += " temp-parse-"+str(counter)+".txt"
   threads[counter].start() 
   #del(threads[counter].lines)

#del(a)
#gc.collect()
   
#f.close()
cat += " >> parsed_"+filename.split('/')[-1]

all = []

for x in range(1, counter+1):
   threads[x].join()
   print x, "joined"
os.system(cat)
os.system("rm -f temp-parse-*")  

fin = open("parsed_"+filename.split('/')[-1])
fin.readline()
bases = 0
 
for l in fin:
   bases+=1

fin.close()

print "parsed_"+filename.split('/')[-1]
print "Bases:\t"+str(bases)
now = time.time()-start
print int(now/60), int(now%60)

