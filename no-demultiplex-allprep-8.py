#! /usr/bin/env python2.6

import os, sys, math, time, operator
import subprocess
from itertools import imap
from optparse import OptionParser
from collections import defaultdict

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2015
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------
#
#
#For a quick view of parameters from the command line, simply use the -h for help.
#no-demultiplex-allprep-8.py -h
#
#This program that looks at all .fastq or .fq files in the current working directory (lib#.fq)
# files, does optional combinations of: 'N' filtering, primary and
#secondary adapter contamination, quality conversion, quality trimming (mean quality of 20 over 5 base window),
#length trimming (default 35)
#
#Input:
#This script works on a directory of fastq files, the files must end with a common fastq file extension
#For each file, depending on the parameters will process a predemultiplexed read file

#Default Behavior
#Rejected reads not output
#'N' nuecleotide reads removed
#Adapter contamination trimmed out
#Qualtites are Phred+33 (samger or illumina 1.8+)
#minimum trimmed sequence length is 35
#single ended reads
#prefix for copy of file that have been changed will be "prepped"

#The options to change behavior
#-E, output rejected reads per file
#-n, allow reads with 'N' nuecleotides
#-a, do not trim out adapter conamination
#-q, convert form phred+64 (illumina 1.5-1.7) to phred+33(samger, 1.8+)
#-M <n>, change minimum sequence length to integer <n>
#-p, ALL input libs are pair ended. PE filss must be pre-interleaved
#-o, a prefix to apply to the new processed file of reads, this work will be put on the front of all of the result files


usage = "\npath/%prog -h"
parser = OptionParser(usage=usage)
parser.add_option("-E", "--error", dest="errfile", action="store_true", default = False, help="Create a rejected reads file.")
parser.add_option("-n", "--allowN", dest="ncheck", action="store_true", default = False, help="Do not reject any sequence with'N' nucleotide automatically.")
parser.add_option("-a", "--allowAdapter", dest="acheck", action="store_true", default = False, help="Do not reject any sequence with adapter contamination automatically.")
parser.add_option("-q", "--qualities", dest="qual", action="store_true", default = False, help="Convert Illumina 1.5 read qualities to sanger.")
parser.add_option("-M", "--minseqlen", dest="minSeqLen", type = "int", default=35, help="Minimum sequence Length")
parser.add_option("-p", "--paired", dest="paired", action="store_true", default = False, help="Create a rejected reads file.")
parser.add_option("-o", "--outprefix", dest='o', default="prepped", help="Output prefix")
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



#used to output to rejected file
def error_out(efile, read, error):
   if opt.errfile != False:
      eout = read[0][0]+read[0][1]+"+"+error+'\n'+read[0][3]
      if opt.paired != False:
         eout += ''.join(read[1])
      efile.write(eout)


#Adapter to look for
#secondadapt = "AGATCGGAAG"
#mainadapt = "AGATCGGAAG"


li = os.listdir(os.getcwd())
reduced = filter(lambda x: (x.endswith(".fq") or x.endswith(".fastq") or x.endswith(".FQ") or x.endswith(".FASTQ")) and "rejected-" not in x and opt.o not in x, li)

mainadapt = "AGATCGGAAGAGC"


def get_len(filen):
   #open all command parameter indicated read files   
   try:
      num = 4
      if opt.paired == True:
         num = 8
      return(file_len(filen)/num)
   except:
      parser.error("Please specify reads files. Please check your command line paramters with -h or --help")


# go through all attached files one read at a time   

start = time.time()


for one_file in reduced:
   print one_file
   frlen = get_len(one_file)
   ct = 0
   good = 0
   bad = 0
   f = open(one_file)
   output_file = open(opt.o+"-"+one_file, 'w')
   if opt.errfile != False:
      efile = open('rejected-'+opt.o+"-"+one_file, 'w')
   while 1:
      if ct % 20000 == 8:
         status(ct, frlen, good, bad)
      ct+=1
      bad_flag = 0
      error = ''

      #pull reads frm files
      if opt.paired == True:
         forwardname =  f.readline()
         if forwardname == "":
            break
         forwardseq = f.readline()
         p1 = f.readline()
         forwardqual = f.readline()
         if opt.qual == True:
            forwardqual = "".join(map((lambda x: chr(ord(x)-31)),forwardqual[:-1]))+'\n'
         reversename =  f.readline()
         reverseseq = f.readline()
         p2 = f.readline()
         reversequal = f.readline()
         if opt.qual == True:
            reversequal = "".join(map((lambda x: chr(ord(x)-31)),reversequal[:-1]))+'\n'
         read = [[forwardname, forwardseq, p1, forwardqual], [reversename, reverseseq, p2, reversequal]]
         if '' in read[1]:
            x[3323]
      else:
         forwardname =  f.readline()
         if forwardname == "":
            break
         forwardseq = f.readline()
         p1 = f.readline()
         forwardqual = f.readline()
         if opt.qual == True:
            forwardqual = "".join(map((lambda x: chr(ord(x)-31)),forwardqual[:-1]))+'\n'
         read = [[forwardname, forwardseq, p1, forwardqual], []]         

      #check for 'N's in sequence
      if opt.ncheck == False:
         if 'N' in forwardseq:
            bad_flag = 1 
            error += "NinF"
         if opt.paired != False:
            if 'N' in reverseseq:
               bad_flag = 1 
               error += "NinR"
         if bad_flag != 0:
            if opt.errfile != False:
               error_out(efile, read, error)
            continue
   
               
      #check for main type of adapter contamination and trim as needed
      if opt.acheck == False:   
         if mainadapt in forwardseq:
            t1 = forwardseq[:forwardseq.index(mainadapt)]
            l1 = len(t1)          
            forwardseq = t1+'\n'
            forwardqual = forwardqual[:l1]+'\n'
         if opt.paired != False:
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
         if opt.paired != False:       
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
      if opt.paired != False:
         quals2 = map(lambda x: ord(x)-33, reversequal[:-1])
         for x in range(len(quals2)-4):
            cut2 = quals2[x:x+5]
            ave2 = float(sum(cut2))/5.0
            if ave2 < minMeanQual:
               reversequal = reversequal[:x]+'\n'
               reverseseq = reverseseq[:x]+'\n'
               break
   
      #check that chopped length is not too small
      if len(forwardseq) < opt.minSeqLen+1 or len(forwardqual) < opt.minSeqLen+1:
         bad_flag = 1
         error += "F too short"
      if opt.paired != False: 
         if len(reverseseq)< opt.minSeqLen+1 or len(reversequal) < opt.minSeqLen+1:
            bad_flag = 1
            error += "R too short" 
            x[3333]
      if bad_flag != 0:
         bad+=1
         if opt.errfile != False:
            error_out(efile, read, error)
         continue
   
      outline = forwardname+forwardseq+'+\n'+forwardqual
      if opt.paired != False:
         outline += reversename+reverseseq+'+\n'+reversequal
         
      output_file.write(outline)
      good+=1
   status(ct, frlen, good, bad)   
   print ""   
   if opt.errfile != False:
     efile.close()
   
      
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
   
   
   

      











