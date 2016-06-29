#! /usr/bin/env python2.6
import os, sys, math, random
from optparse import OptionParser
from collections import defaultdict

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-f", "--file", dest="f", help="Input MAPS2 mutation file.")
parser.add_option("-r", "--ref", dest="r", help="Input reference file.")
parser.add_option("-o", "--outfile", dest="o", help="output file")
(opt, args) = parser.parse_args()

random.seed(8)

def readFasta(filename):
   f = open(filename)       
   #read in a fasta
   name = ''
   seq = ''
   seqs = {} 
   while True:
      l = f.readline()
      if l == '':
         seqs[name] = seq
         break
      elif l[0] == '>':
         if name != '':
            seqs[name] = seq

         seq = ''
         name = l[1:-1]
      else:
         seq+=str.upper(l).replace('\n','')
   f.close()
   
#KVK   print seqs.keys()

   return seqs



f = open(opt.f)   
ref = readFasta(opt.r)  
frame = {}

#KVK print seqs['IWGSC_CSS_5BL_scaff_10909027']

frame['hom'] = {}
frame['het'] = {}

f.seek(0)
f.readline()
for l in f:
   x = l[:-1].split('\t')
   chrom = x[0]
   pos = int(x[1])
   refb = x[4]

#KVK   print chrom, "\t", pos, "\n"

   if '*' in refb or '+' in refb:
      continue
# deleted rice-specific samples 
   if x[6] in ['Kronos0']:
      print "Kronos0 found", x
      continue
   ctype = x[10]
   hohe = x[7]   

   tseta = ref[chrom][pos-51:pos-41]
   tsetaIndex = [i for i, y in enumerate(tseta) if y == refb]
   if len(tsetaIndex) == 0:
      seta = ''.join(map(lambda y: '.', range(41))) 
   else:
      aRefIndex = random.choice(tsetaIndex)
      randSetA = ref[chrom][pos-51+aRefIndex-20:pos-51+aRefIndex+21]
      seta = randSetA
   tsetc = ref[chrom][pos+40:pos+49]
   tsetcIndex = [i for i, y in enumerate(tsetc) if y == refb]
   if len(tsetcIndex) == 0:
      setc = ''.join(map(lambda y: '.', range(41)))
   else:
      cRefIndex = random.choice(tsetcIndex)
      randSetC = ref[chrom][pos+40+cRefIndex-20:pos+40+cRefIndex+21]
      setc = randSetC

   setb = ref[chrom][pos-21:pos-1]+refb+ref[chrom][pos:pos+20]
   allset = seta + setb + setc
  
   if len(allset) != 123:	# Paul B - 11.5.2016 - allset looks to be all the positions for the analysis. They should add up to 123. If they don't then I think the row should ignored!  
      print x				# these will be the rows that can't provide 123 positions because they are at the start or end of the contig!
      #Paul B. -11.5.2016 - added "continue" to test the Cadenza D genome because it kept crashing  
      continue
#KVK commented out line below      
#     raise Exception("Border Case: main array size too short. Missing values on sides")
   if 'N' in allset:
      continue
   if 'S' in allset:
      continue
   if 'R' in allset:
      continue   
   if 'Y' in allset:
      continue
   if 'W' in allset:
      continue
   if 'M' in allset:
      continue
   if 'K' in allset:
      continue
   # Skip other ambiguity codes in the reference otherwise script crashes (Paul B. added 10.2.2016): 
   if 'V' in allset or 'D' in allset or 'B' in allset or 'H' in allset:
      continue    
   if ctype not in frame[hohe]:
      frame[hohe][ctype] = map(lambda y: {'A': 0, 'C': 0, 'T': 0, 'G': 0}, range(len(allset)))  
   for i in range(len(allset)):
      base = allset[i]
      if base != '.':
         frame[hohe][ctype][i][base] +=1     
  
f.close()
o = open(opt.o, 'w')
headline = ["Hom/Het", "Type", 'Base']
guide = range(-62,-21)+range(-20,21)+range(21,62)
guide = map(lambda x: str(x), guide)

for i in range(len(allset)):
   ti = guide[i]
   if ti == '0':
      ti = "Ref"
   headline += [ti]

o.write('\t'.join(headline)+'\n')
for thohe in ['hom', 'het']:
   ind1 = frame[thohe].keys()
   ind1.sort()
   for ct in ind1:
      tset = frame[thohe][ct]
      aline = [thohe, ct, 'A']
      cline = [thohe, ct, 'C']
      gline = [thohe, ct, 'G']
      tline = [thohe, ct, 'T']
      for i in range(len(allset)):
         aline.append(str(tset[i]['A']))
         cline.append(str(tset[i]['C']))
         gline.append(str(tset[i]['G']))
         tline.append(str(tset[i]['T']))
      o.write('\t'.join(aline)+'\n')
      o.write('\t'.join(cline)+'\n')
      o.write('\t'.join(gline)+'\n')
      o.write('\t'.join(tline)+'\n')

o.close()
   
   
   
   
   
   
  
   
   
