#!/usr/bin/env python
#fix_qv -- fix qvalue after combining LSA and LAA results

#License: BSD

#Copyright (c) 2008 Li Charles Xia
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions
#are met:
#1. Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#3. The name of the author may not be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
#IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
#IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
#INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
#NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
#THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#public libs
import sys, csv, re, os, time, argparse, string, tempfile, subprocess
#numeric libs
import numpy as np
import scipy as sp
try:
  #debug import
  from . import lsalib
except ImportError:
  #install import
  from lsa import lsalib
  #np.seterr(all='raise')

#assume column index Pi and Qi is given

#input arguments:
#raw_file, fixed_file, Pi, Qi

ws=os.path.join(os.environ.get("HOME"),'tmp','multi')
print("tmpDir=",ws)

def main():  

  parser = argparse.ArgumentParser(description="Fix Qvalue Commandline tool for LSA and LA")

  parser.add_argument("rawInput", metavar="rawInput", type=argparse.FileType('r'), help="the raw input file")
  parser.add_argument("fixedOutput", metavar="fixedOutput", type=argparse.FileType('w'), help="the fixed output file")
  parser.add_argument("-pi", "--pIndex", dest="pIndex", type=int, default=7, help="pvalue column index")
  parser.add_argument("-qi", "--qIndex", dest="qIndex", type=int, default=8, help="qvalue column index")

  arg_namespace=parser.parse_args()
  rawInput=vars(arg_namespace)['rawInput']
  fixedOutput=vars(arg_namespace)['fixedOutput']
  qIndex=vars(arg_namespace)['qIndex']
  pIndex=vars(arg_namespace)['pIndex']

  #read contents
  full_table = []
  for row in csv.reader(rawInput, delimiter="\t"):
    full_table.append(row)
  
  #collect pvalues
  pvalues = []
  j = 0
  for row in full_table:
    if j ==0:
      j=j+1
      continue
    else:
      #print len(row)
      pvalues.append(row[pIndex-1])

  pvalues=np.array(pvalues, dtype='float')
  #print pvalues
  
  #pvalues to qvalues
  qvalues = lsalib.storeyQvalue(pvalues)
  
  #insert qvalues
  j=0
  for row in full_table:
    if j ==0:
      j=j+1
      continue
    else:
      #print len(row)
      row[qIndex-1]=qvalues[j-1]
      j=j+1
  
  #write contents
  write_table = csv.writer(fixedOutput,delimiter="\t")
  write_table.writerows(full_table) 

if __name__ == "__main__":
  main()
