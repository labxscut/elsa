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
import sys, csv, re, os, time, argparse, string, tempfile, subprocess
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

ws=os.path.join(os.environ.get("HOME"),'tmp','multi')
print "tmpDir=",ws

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

def gen_singles(multiInput, multiOutput):
  header, content=get_content(multiInput)
  multiname=multiInput.name
  i=1
  singles=[]
  results=[]
  ends=[]
  for line in content:
    singlename=multiname+".%d" % i
    tmp=open(os.path.join(ws, singlename), "w")
    print >>tmp, header.rstrip('\n')
    print >>tmp, line.rstrip('\n')
    tmp.close()
    singles.append(tmp.name) 
    results.append(tmp.name+".tmp") 
    ends.append(tmp.name+".end")
    i+=1
  return singles, results, ends

PBS_PREAMBLE = """#!/bin/bash
#PBS -N %s 
#PBS -S /bin/bash 
#PBS -j oe
#PBS -o %s
#PBS -l walltime=299:00:00
#PBS -l nodes=1:ppn=1,mem=%d000mb,vmem=%d000mb,pmem=%d000mb"""
vmem=2
print "vmem=", str(vmem)+"000mb"

def gen_pbs(singleFile, singleCmd, workDir, singleEnd):
  singleResult=singleFile+".tmp"
  singlePBS=open(singleFile+".pbs", "w")
  print >>singlePBS, PBS_PREAMBLE % (os.path.basename(singleFile), os.path.basename(singleFile)+".log", vmem, vmem, vmem)
  print >>singlePBS, "cd %s" % workDir
  print >>singlePBS, singleCmd % (singleFile, singleResult)
  print >>singlePBS, "touch %s" % singleEnd
  singlePBS.close()
  return singlePBS.name

def gen_output(multiOutput, resultFiles):
  i=0
  for resultFile in resultFiles:
    header, content = get_content(resultFile)
    if i==0:
      print >>multiOutput, "".join([header]+content)
      i+=1
    else:
      print >>multiOutput, "".join(content)
  return

def ssa_pbs(singlePBS):
  try:
    tmp=subprocess.Popen("ssa.py %s" % singlePBS, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    print "submitted", singlePBS
  except ValueError:
    quit()
  return tmp[0]

def main():  

  parser = argparse.ArgumentParser(description="Multiline Input Split and Combine Tool for LSA and LA")

  parser.add_argument("multiInput", metavar="multiInput", type=argparse.FileType('r'), help="the multiline input file")
  parser.add_argument("multiOutput", metavar="multiOutput", type=argparse.FileType('w'), help="the multiline output file")
  parser.add_argument("singleCmd", metavar="singleCmd", help="single line command line in quotes")
  parser.add_argument("workDir", metavar="workDir", help="set current working directory")

#  """la_compute ARISA.depCmax.txt ARISA.depCmax.S5_L75_Ptheo.lsaq ARISA.depCmax.S5_L75_Ptheo.la -s 114 -r 1 -p 1000"""
  arg_namespace=parser.parse_args()
  multiInput=vars(arg_namespace)['multiInput']
  multiOutput=vars(arg_namespace)['multiOutput']
  singleCmd=vars(arg_namespace)['singleCmd']
  workDir=vars(arg_namespace)['workDir']

  singleFiles,resultFiles,endFiles=gen_singles(multiInput,multiOutput)
  for endFile in endFiles:
    if os.path.exists(endFile):
      os.remove(endFile)

  inProgress=set()
  endJob=set()
  while(len(singleFiles)!=0):
    singleFile=singleFiles.pop()
    endFile=endFiles.pop()
    pbsFile=gen_pbs(singleFile, singleCmd, workDir, endFile)
    inProgress.add(endFile)
    ssa_pbs(pbsFile)

  while(len(endJob)!=len(inProgress)):
    for job in inProgress:
      time.sleep(10)
      if os.path.exists(job):
        header, content = get_content(open(job[:-4]+".tmp",'r'))
        if len(endJob)==0:
          print >>multiOutput, "".join([header]+content)
        else:
          print >>multiOutput, "".join(content)
        os.remove(job)
        endJob.add(job)
        print "ended", job
        print "remaining jobs", inProgress.difference(endJob), "total", len(inProgress.difference(endJob))

  #gen_output(multiOutput, resultFiles)

if __name__ == "__main__":
  main()
