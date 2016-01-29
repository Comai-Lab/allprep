#! /usr/bin/env python

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, 2011
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------

import os, sys, math
from optparse import OptionParser

#interleaveSwitcher.py

#Usage: programName.py i forwardFile.fq reverseFile.fq outFileName.fq
#       ______OR_______
#Usage: programName.py u interleavedFile

#This program has two operating modes, defined by the first command line parameter.

#Interleave mode, first parameter an 'i' or an 'I'.
#   This will take two read files of equal length and interleave them.
#Input Parameters:
#forwardFile.fq = original file input, sanger or illumina file forward side
#reverseFile = original file input, sanger or illumina reverse side
#outFileName.fq = interleaved output file name/location

#Uninterleave mode, first parameter an 'u' or an 'U'.
#   This will take one interleaved paired read file and split them into two paired read files.
#   The -1.fq and -2.fq result files will be forward and backward respectively.
#Input Parameters:
#interleavedFile.fq = original file input, sanger or illumina interleaved paired reads

#Be sure to put the i/u as the first parameter!

#--------------------------------------------------------------------------

usage = "\npath/%prog -h"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file1", dest="f", help="Forward read .fq file OR already interleaved file")
parser.add_option("-r", "--reverse", dest="r", default = False, help="Reverse read .fq file if interleaving")
parser.add_option("-o", "--output", dest="o", default = False, help="output file name for interleaved OR prefix for two output uninterleaved")
(opt, args) = parser.parse_args()

#get operating mode
if opt.r != False:
   mode = 'i'
else:
   mode = 'u'

#if interleaving
if mode == 'i':
   f = open(opt.f)
   r = open(opt.r)
   o = open(opt.o, 'w')
   
   #read in one line from both files at a time and output together
   while True:
      name1 = f.readline()
      if name1 == "":
         break
      seq1 = f.readline()
      p1 = f.readline()
      qual1 = f.readline()
      name2 = r.readline()
      seq2 = r.readline()
      p2 = r.readline()
      qual2 = r.readline()
      o.write(name1+seq1+p1+qual1+name2+seq2+p2+qual2)
   
   f.close()
   r.close()
   o.close()
   
   
#if uninterleaving
#create two files for the paired data
#then read in 8 lines and sent the first 4 to forward and 2nd 4 to reverse
#assuming fastq or illumina 4 line read format
elif mode == 'u':
   u = open(sys.argv[2])
   name = opt.o
   f = open(name+"-1.fq",'w')
   r = open(name+"-2.fq",'w')

   while True:
      name1 = u.readline()
      if name1 == '':
         break
      seq1 = u.readline()
      p1 = u.readline()
      qual1 = u.readline()
      name2 = u.readline()
      seq2 = u.readline()
      p2 = u.readline()
      qual2 = u.readline()
      f.write(name1+seq1+p1+qual1)
      r.write(name2+seq2+p2+qual2)
   
   f.close()
   r.close()
   u.close()
#if invalid mode specified
else:
   print "Use correct run format, see script header for details"



