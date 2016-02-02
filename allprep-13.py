#! /usr/bin/env python2.6

import os, sys, math, time, operator
import subprocess
from itertools import imap
from optparse import OptionParser
from collections import defaultdict

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2016
#This work is the property of UC Davis Genome Center - Comai Lab
#This is shared under a Creative Commons BY-NC-ND 4.0 license
#https://creativecommons.org/licenses/by-nc-nd/4.0/

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------
#
#Usage: 
#For all modes the run command looks like this, with [...] indicating files needed by specific read type, and {....} indication optional parameters
#allprep-8.py -b barcode-file.txt -f forward-read-file.fq [-r reverse-read-file.fq] [-i index1-file.fq] [-I index2-file.fq] {-m} {-E} {-n} {-q} {-d}
#
#For a quick view of parameters form the command line, simply use the -h for help.
#allprep-8.py -h
#
#This program that takes a barcode file and splits the lane sequence.txt files into specified
#library.txt (lib#.txt) files, does barcode check and match, 'N' filtering, primary and
#secondary adapter contamination, quality conversion, quality trimming (mean quality of 20 over 5 base window),
#length trimming (default 35), and library separation.
#
#Input:
#This script takes a barcode file as specified in the sample sheet in the README, as well as 
#the forward read file, and optionally the reverse, index and secondary index
#These files are loaded with -b, -f, -r, -i, and -I respectively.
#There are also five additional operating options:
#-m for mismatch mode, this allows a 1 bp mismatch between the index / barcode and the best matching barcode. 
#-E for error reads mode, outputs the rejected reads to a file
#-n for N allowed mode, this turns off the check that rejects a read if there are any 'N' nucleotides in the sequence
#-q to convert a file using Illumina 1.5 qualities to Sanger/Illumina1.8 (standard)
#-D for only demultiplexing mode, this only splits the reads by barcode and does no additional trimming or checks
#
#This program can take a file with most combinations of barcode/indexing. The possibilities are:
#1. single ended, barcoded (read type 0)
#2. pair ended, barcoded (read type 1)
#3. single ended, one index (read type 2)
#4. single ended, two index (read type 5)
#5. pair ended, one index (read type 3)
#6. pair ended, two index (read type 4)
#
#Please see the README-barcode-file.txt and sample-barcode-file.txt for help in creating the barcode file for your dataset.


usage = "\npath/%prog -h"
parser = OptionParser(usage=usage)

parser.add_option("-b", "--barcode", dest="bfile", default = False, help="Barcode/Index file with all lib information. See sample table for details")
parser.add_option("-f", "--forward", dest="f", help="Forward read .fq file")
parser.add_option("-r", "--reverse", dest="r", default = False, help="Reverse read .fq file")
parser.add_option("-i", "--index", dest="i", default = False, help="Index for forward read, or index for paired read if only one index")
parser.add_option("-I", "--indexreverse", dest="i2", default = False, help="Reverse read index file.")
parser.add_option("-m", "--mismatch", dest="miss", action="store_true", default = False, help="Allow one base difference for indexed barcodes")
parser.add_option("-E", "--error", dest="errfile", action="store_true", default = False, help="Create a rejected reads file.")
parser.add_option("-n", "--allowN", dest="ncheck", action="store_true", default = False, help="Do not reject any sequence with'N' nucleotide automatically.")
parser.add_option("-N", "--Ncut", dest="ntrim", action="store_true", default = False, help="Do Trim reads by N instead f rejecting them.")
parser.add_option("-q", "--qualities", dest="qual", action="store_true", default = False, help="Convert Illumina 1.5 read qualities to sanger.")
parser.add_option("-D", "--demultiplex", dest="demulti", action="store_true", default = False, help="ONLY Demultiplex, not trimming or checking at all")
parser.add_option("-M", "--minseqlen", dest="minSeqLen", type = "int", default=35, help="Minimum sequence Length")


(opt, args) = parser.parse_args()

minMeanQual = 20.0



# Uses wc to get te number of lines in the file
def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

#prints out a status bar
def status(cur, total, good, bad):
      now = time.time()-start
      perdone = cur/float(total-1)
      percent = (perdone)*100
      emin = now/60
      esec = now%60
      rmin = (total-cur)/(cur/now)/60
      rsec = (total-cur)/(cur/now)%60
      sys.stdout.write('\r')      
      sys.stdout.write(">%-20s< %i%% Elapsed: %im %is Remaining: %im %is Good: %i Rejected %i Total: %i   " % ('='*int((perdone*20)),percent,emin,esec,rmin,rsec, good, bad, good+bad))
      sys.stdout.flush()

#complement sequences
def comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'R': 'R', 'Y': 'Y', 'S': 'S', 'W': 'W', 'M': 'M', 'K': 'K'}
    complseq = [complement[base] for base in seq]
    return complseq

#reverse complemtent
def rev(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(comp(seq))

#gets 4 reads from all given command parameter files
def get_reads(readset):
   ind = 0
   reads = [[],[],[],[]]
   for fh in readset:
      if fh != '':
         name = fh.readline()
         if name == '':
            break
         seq = fh.readline()
         plus = fh.readline()
         if opt.qual == True:
            qualt = fh.readline()
            qual = "".join(map((lambda x: chr(ord(x)-31)),qualt[:-1]))+'\n'
         else:
            qual = fh.readline()
         reads[ind] = [name, seq, plus, qual]
      ind+=1
   if reads != [[],[],[],[]]:
      return reads
   else:
      return "EOF"

#extracts sequences from all possible reads
def get_codes(readset):
   try:
      code1 = read[0][1][:-1]
   except:
      code1 = ''
   try:
      code2 = read[1][1][:-1]
   except:
      code2 = ''
   try:
      codei1 = read[2][1][:-1]
   except:
      codei1 = ''
   try:
      codei2 = read[3][1][:-1]
   except:
      codei2 = '' 
   if lenindex != 0 and opt.i == False:
      codei1 = read[0][0].split('#')[-1].split('/')[0]
      try:
         codei2 = read[1][0].split('#')[-1].split('/')[0]
      except:
         codei2 = ''
   return [code1, code2, codei1, codei2]  

#checks start of sequences against possible barcodes, returns a file handle, the code type, and the barcode
def check_codes(curcodes):
   cfh = ''
   mtype = ''
   if len(codes['pindex2']) != 0:
      for poss in codes['pindex2'].keys():
         pcode1, pcode2 = poss.split('-')
         if opt.miss == False:
            if pcode1 == curcodes[2][:len(pcode1)] and pcode2 == curcodes[3][:len(pcode2)]:
               cfh = codes['pindex2'][poss]
               mtype = 'pindex2'
               break
         else:
            d1 = sum(imap(ne, pcode1, curcodes[2]))
            d2 = sum(imap(ne, pcode2, curcodes[3]))
            if d1 <= 1 and d2 <= 1:
               cfh = codes['pindex2'][poss]
               mtype = 'pindex2'
               break
   if len(codes['sindex2']) != 0 and cfh == '':
      for poss in codes['sindex2'].keys():
         pcode1, pcode2 = poss.split('-')
         if opt.miss == False:
            if pcode1 == curcodes[2][:len(pcode1)] and pcode2 == curcodes[3][:len(pcode2)]:
               cfh = codes['sindex2'][poss]
               mtype = 'sindex2'
               break
         else:
            d1 = sum(imap(ne, pcode1, curcodes[2]))
            d2 = sum(imap(ne, pcode2, curcodes[3]))
            if d1 <= 1 and d2 <= 1:
               cfh = codes['sindex2'][poss]
               mtype = 'sindex2'
               break
   if len(codes['pindex1']) != 0 and cfh == '':
      for poss in codes['pindex1'].keys():
         if opt.miss == False:
            if poss == curcodes[2][:len(poss)]:
               cfh = codes['pindex1'][poss]
               mtype = 'pindex1'
               break
         else:
            d1 = sum(imap(ne, poss, curcodes[2]))  
            if d1 <= 1:
               cfh = codes['pindex1'][poss]
               mtype = 'pindex1'
               break  
   if len(codes['sindex1']) != 0 and cfh == '':
      for poss in codes['sindex1'].keys():
         if opt.miss == False:
            if poss == curcodes[2][:len(poss)]:
               cfh = codes['sindex1'][poss]
               mtype = 'sindex1'
               break
         else:
            d1 = sum(imap(ne, poss, curcodes[2]))  
            if d1 <= 1:
               cfh = codes['sindex1'][poss]
               mtype = 'sindex1'
               break
   if len(codes['old']) != 0 and cfh == '':
      for poss in codes['old'].keys():
         tposs = rev(poss)+overhangs[poss]
         if opt.miss == False:
            if tposs == curcodes[0][:len(tposs)] and (tposs == curcodes[1][:len(tposs)] or opt.r == False):
               cfh = codes['old'][poss]
               mtype = 'old'
               break
         else:
            d1 = sum(imap(ne, tposs, curcodes[0]))
            d2 = sum(imap(ne, tposs, curcodes[1]))
            if d1 <= 1 and (d2 <= 1 or opt.r == False):
               cfh = codes['old'][poss]
               mtype = 'old'
               break
   return [cfh, mtype, poss]  

#used to output to rejected file
def error_out(read, error):
   if opt.errfile != False:
      eout = read[0][0]+read[0][1]+"+"+error+'\n'+read[0][3]
      if opt.r != False:
         eout += ''.join(read[1])
      if opt.i != False:
         eout += ''.join(read[2])
      if opt.i2 != False:
         eout += ''.join(read[3])
      efile.write(eout)


#Adapter to look for
#secondadapt = "AGATCGGAAG"
#mainadapt = "AGATCGGAAG"



mainadapt = "AGATCGGAAGAGC"


ne = operator.ne
######### READ IN BARCODE FILE
try:
   bt = open(opt.bfile)
   bt.readline()
except:
   parser.error("Missing valid barcode table. Please check your command line paramters with -h or --help")
   
codes = {}
codes['old'] = {}
codes['sindex1'] = {}
codes['sindex2'] = {}
codes['pindex1'] = {}
codes['pindex2'] = {}
rembases = {}
overhangs = {}

handles={}
#read in the barcode table, sorting by type
for l in bt:
  
   x= l[:-1].split('\t')
   if x == ['']:
      continue
   x = map(lambda y: y.replace(' ',''), x)
   libname, rtype, code1, code2, over, rembase = x
   rembase = rembase.replace('.','0')
   if int(rtype) == 0 or int(rtype) == 1:
      ctype = 'old'
      name = code1
   if int(rtype) == 2:   
      ctype = 'sindex1'
      name = code1
   if int(rtype) == 3:
      ctype = 'pindex1'
      name = code1
   if int(rtype) == 4:
      ctype = 'pindex2'
      name = code1+'-'+code2
   if int(rtype) == 5:
      ctype = 'sindex2'
      name = code1+'-'+code2      
   codes[ctype][name] = libname
   if libname not in handles:
      handles[libname] = open(libname+'.fq', 'w')
   else:
      print libname
      parser.error("Duplicate Filename !!! :"+libname)
   rembases[name] = int(rembase)
   over = over.replace('.','')
   overhangs[name] = over

bt.close()
   
lenindex = len(codes['sindex1']) + len(codes['pindex1']) + len(codes['pindex2']) + len(codes['sindex2'])

#duh, do not run for both  
if opt.miss == True:
   if len(codes['old']) != 0 and lenindex != 0:
      parser.error("Missmatch mode cannot be run with both barcodes and indexed reads. Please check your command line paramters with -h or --help")
      
#open all command parameter indicated read files   
try:
   frlen = file_len(opt.f)/4
   #frlen = 296146516
   ffile = open(opt.f)
except:
   parser.error("Please specify reads files. Please check your command line paramters with -h or --help")
if opt.r !=  False:
   rfile = open(opt.r)
else:
   rfile = ''
if opt.i !=  False:
   ifile = open(opt.i)
else:
   ifile = ''
if opt.i2 !=  False:
   i2file = open(opt.i2)
else:
   i2file = ''
if opt.errfile != False:
   efile = open("rejected-reads-"+opt.bfile, 'w')
creadset = [ffile, rfile, ifile, i2file]   

# go through all attached files one read at a time   
ct = 0
good = 0
bad = 0
start = time.time()
while 1:
   if ct % 20000 == 8:
      status(ct, frlen, good, bad)
   ct+=1

   #pull reads frm files
   read = get_reads(creadset)
   if read == "EOF":
      break
      
   #get relevant barcodes/indexes
   curcodes = get_codes(read)
   #valid_codes = check_codes(codes)
   
   cfh, mtype, poss = check_codes(curcodes)
   bad_flag = 0
   error = ''
   if cfh == '':
      bad_flag = 1
      error = "No Barcode Match"
      bad+=1
      error_out(read, error)
      continue
   
   nameflag = 0
   if opt.i == False:
      nameflag = 1

   if mtype == 'old':
      forwardname = read[0][0][:-1]+'#'+poss+'\n'
      forwardseq = read[0][1][len(poss)+len(overhangs[poss]):]
      forwardqual = read[0][3][len(poss)+len(overhangs[poss]):]
      if opt.r != False:
         reversename = read[1][0][:-1]+'#'+poss+'\n'
         reverseseq = read[1][1][len(poss)+len(overhangs[poss]):]
         reversequal = read[1][3][len(poss)+len(overhangs[poss]):]
   elif mtype == 'sindex1' or mtype == 'pindex1':
      if nameflag == 0:
         forwardname = read[0][0][:-1]+poss+'\n'
      else:
         forwardname = read[0][0][:-1]+'\n'
      forwardseq = read[0][1]
      forwardqual = read[0][3]
      if mtype == 'pindex1':
         if nameflag == 0:
            reversename = read[1][0][:-1]+poss+'\n'
         else:
            reversename = read[1][0][:-1]+'\n'
         reverseseq = read[1][1]
         reversequal = read[1][3]
   elif mtype == 'pindex2':  
      poss1, poss2 = poss.split('-')     
      if nameflag == 0:
         forwardname = read[0][0][:-1]+poss1+'\n'
         reversename = read[1][0][:-1]+poss2+'\n'
      else:
         forwardname = read[0][0][:-1]+'\n'
         reversename = read[1][0][:-1]+'\n'
      forwardseq = read[0][1]
      forwardqual = read[0][3]

      reverseseq = read[1][1]
      reversequal = read[1][3]
   elif mtype == 'sindex2':  
      poss1, poss2 = poss.split('-')    
      if nameflag == 0: 
         forwardname = read[0][0][:-1]+poss1+':'+poss2+'\n'
      else:
         forwardname = read[0][0][:-1]+'\n'
      forwardseq = read[0][1]
      forwardqual = read[0][3]
   else:
      parser.error("Error: Check barcode splitting")
   
   #remove bases
   if rembases[poss] != 0:
      forwardseq = forwardseq[rembases[poss]:]
      forwardqual = forwardqual[rembases[poss]:]
      if opt.r != False:
         reverseseq = reverseseq[rembases[poss]:]
         reversequal = reversequal[rembases[poss]:]

   if opt.demulti == False:
             
      #check for main type of adapter contamination and trim as needed
   
      if mainadapt in forwardseq:
         t1 = forwardseq[:forwardseq.index(mainadapt)]
         l1 = len(t1)          
         forwardseq = t1+'\n'
         forwardqual = forwardqual[:l1]+'\n'
      if opt.r != False:
         if mainadapt in reverseseq:
            t1 = reverseseq[:reverseseq.index(mainadapt)]
            l1 = len(t1)          
            reverseseq = t1+'\n'
            reversequal = reversequal[:l1]+'\n'       
         
      if rev(mainadapt) in forwardseq:
         t1 = forwardseq[:forwardseq.index(rev(mainadapt))]
         l1 = len(t1)          
         forwardseq = t1+'\n'
         forwardqual = forwardqual[:l1]+'\n'  
      if opt.r != False:       
         if rev(mainadapt) in reverseseq:
            t2 = reverseseq[:reverseseq.index(rev(mainadapt))]
            l2 = len(t2)          
            reverseseq = t2+'\n'
            reversequal = reversequal[:l2]+'\n'     
   
    
      #calculate avg quality and trim if lower then threshold (20)
   
      quals = map(lambda x: ord(x)-33, forwardqual[:-1])
      for x in range(len(quals)-4):
         cut = quals[x:x+5]
         ave = float(sum(cut))/5.0
         if ave < minMeanQual:
            forwardqual = forwardqual[:x]+'\n'
            forwardseq = forwardseq[:x]+'\n'
            break
      if opt.r != False:
         quals2 = map(lambda x: ord(x)-33, reversequal[:-1])
         for x in range(len(quals2)-4):
            cut2 = quals2[x:x+5]
            ave2 = float(sum(cut2))/5.0
            if ave2 < minMeanQual:
               reversequal = reversequal[:x]+'\n'
               reverseseq = reverseseq[:x]+'\n'
               break

      if opt.ntrim == True:
         if opt.ncheck == False:
            for x in range(len(forwardseq)):
               if forwardseq[x] == 'N':
                  forwardqual = forwardqual[:x]+'\n'
                  forwardseq = forwardseq[:x]+'\n'
                  break
            if opt.r != False:
               for x in range(len(reverseseq)):
                  if reverseseq[x] == 'N':
                     reversequal = reversequal[:x]+'\n'
                     reverseseq = reverseseq[:x]+'\n'
                     break      
      elif opt.ntrim == False:         
         #check for 'N's in sequence
         if opt.ncheck == False:
            if 'N' in forwardseq:
               bad_flag = 1 
               error += "NinF"
            if opt.r != False:
               if 'N' in reverseseq:
                  bad_flag = 1 
                  error += "NinR"
            if bad_flag != 0:
               bad+=1
               error_out(read, error)
               continue
         
   
      #check that chopped length is not too small
      if len(forwardseq) < opt.minSeqLen+1 or len(forwardqual) < opt.minSeqLen+1:
         bad_flag = 1
         error += "F too short"
      if opt.r != False: 
         if len(reverseseq)< opt.minSeqLen+1 or len(reversequal) < opt.minSeqLen+1:
            bad_flag = 1
            error += "R too short" 
      if bad_flag != 0:
         bad+=1
         error_out(read, error)
         continue

   outline = forwardname+forwardseq+'+\n'+forwardqual
   if opt.r != False and mtype != 'sindex1' and mtype != 'sindex2':
      outline += reversename+reverseseq+'+\n'+reversequal
      
   handles[cfh].write(outline)
   good+=1
status(ct, frlen, good, bad)   
print ""   
   
#close all library file handles 
for ctype in handles.keys():
   handles[ctype].close()
         

ffile.close()

if opt.r !=  False:
   rfile.close()
   
if opt.i !=  False:
   ifile.close()
 
if opt.i2 !=  False:
   i2file.close()

if opt.errfile != False:
   efile.close()        
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
   
   
   

      











