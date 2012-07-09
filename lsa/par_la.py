#!/usr/bin/env python
#par_la -- parallel computation of liquid association 

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
import sys, csv, re, os, time, argparse, string, tempfile
#numeric libs
import numpy as np
import scipy as sp
try:
  #debug import
  import lsalib
except ImportError:
  #install import
  from lsa import lsalib
  #np.seterr(all='raise')

#assuming input of form "xxxxx %single_input %single_output xxxxxx", only % to be replaced by single line file
#multiline file format:
#headline
#content1
#content2
#...
#convert to
#file1:
#headline
#content1
#file2:
#headline
#content2
#...
#single result file format
#result_file1
#headline
#result1-1
#result1-2
#...
#result_file2
#headline
#result2-1
#result2-2
#...
#combine back to
#result_file
#headline
#result1-1
#result1-2
#...
#result2-1
#result2-2
#...

#input arguments:
#multiInput, multiOutput, singleCmd

ws=os.path.join( [os.environ.get("HOME"),'tmp','multi'] )

def get_content(file):
  i=0
  header=None
  content=[]
  for line in file:
    if i==0:
      header=line
      i+=1
    else:
      content.append(line)
  return (header, content)

def gen_singles(multiInput):
  header, content=get_content(multiInput)
  multiname=multiInput.name
  i=1
  singles=[]
  results=[]
  for line in content:
    tmp=open(multiname+".%d" % i, "w")
    print >>tmp, header
    print >>tmp, line
    tmp.close()
    singles.append(tmp.name) 
    results.append(tmp.name+".tmp") 
  return singles, results

PBS_PREAMBLE = """#!/bin/bash
#PBS -N %s 
#PBS -S /bin/bash 
#PBS -j oe
##PBS -k oe
#PBS -o %s
#PBS -l walltime=299:00:00
#PBS -l nodes=1:ppn=1,mem=%d000mb,vmem=%d000mb,pmem=%d000mb"""

def gen_pbs(singleFile, singleCmd):
  singleResult=singleFile+".tmp"
  singlePBS=open(singleFile+".pbs", "w")
  print >>singlePBS, PBS_PREAMBLE % (singleFile, singleFile, 1, 1, 1)
  print >>singlePBS, singleCmd % (singleFile, singleResult)
  singlePBS.close()
  return singlePBS.name

def gen_output(multiOutput, resultFiles):
  header = None
  content = []
  for file in resultFiles):
    header, tmp_content = get_content(resultFile)
    content += tmp_content
  print >>multiOutput, "\n".join([header]+content)
  return

def ssa_pbs(singlePBS):
  try:
    print "submit", singlePBS
    #tmp=subprocess.Popen("ssa.py %s" % singlePBS, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
  except OSError:
    quit()
  return

def main():  

  parser = argparse.ArgumentParser(description="Multiline Input Split and Combine Tool for LSA and LA")

  parser.add_argument("multiInput", metavar="multiInput", type=argparse.FileType('r'), help="the multiline input file")
  parser.add_argument("multiOuput", metavar="multiOutput", type=argparse.FileType('r'), help="the multiline output file")
  parser.add_argument("singleCmd", metavar="singleCmd", help="single line command line in quotes")

#  """la_compute ARISA.depCmax.txt ARISA.depCmax.S5_L75_Ptheo.lsaq ARISA.depCmax.S5_L75_Ptheo.la -s 114 -r 1 -p 1000"""
  arg_namespace=parser.parse_args()
  multiInput=vars(arg_namespace)['multiInput']
  multiOutput=vars(arg_namespace)['multiOutput']
  singleCmd=vars(arg_namespace)['singleCmd']

  singleFiles,resultFiles=gen_singles(multiInput)
  while(len(singleFiles)!=0):
    singleFile=singleFiles.pop()
    pbsFile=gen_pbs(singleFile)
    ssa_pbs(pbsFile)

  gen_output(resultFiles)

if __name__ == "__main__":
  main()
