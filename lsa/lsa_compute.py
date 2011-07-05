#!/usr/bin/env python
#lsa-compute -- computation script for LSA package to perform lsa table calculation 

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
  #installed import
  from lsa import lsalib
except ImportError:
  #debug import
  import lsalib
  #np.seterr(all='raise')

def main():  

  # define arguments: delayLimit, fillMethod, permuNum
  parser = argparse.ArgumentParser(description="New LSA Commandline Tool")

  parser.add_argument("dataFile", metavar= "dataFile", type=argparse.FileType('r'), help="the input data file")
  parser.add_argument("resultFile", metavar= "resultFile", type=argparse.FileType('w'), help="the output result file")
  parser.add_argument("-d", "--delayLimit", dest="delayLimit", default=3, type=int, choices=range(0,6),
                    	help="specify the maximum delay possible, default: 3,\n choices: 0 to 6")
  parser.add_argument("-p", "--permuNum", dest="permuNum", default=1000, type=int, choices=[100, 200, 500, 1000, 2000],
                    	help="specify the number of permutations for p-value estimation, default: 1000,\n \
                          choices: 100, 200, 500, 1000, 2000.")
  parser.add_argument("-b", "--bootNum", dest="bootNum", default=100, type=int, choices=[100, 200, 500, 1000, 2000],
                    	help="specify the number of bootstraps for 95%% confidence interval estimation, default: 100,\n \
                          choices: 100, 200, 500, 1000, 2000.")   #use %% to print %
  parser.add_argument("-r", "--repNum", dest="repNum", default=1, type=int,
                    	help="specify the number of replicates each time spot, default: 1,      \n \
                          must be provided and valid. ")
  parser.add_argument("-s", "--spotNum", dest="spotNum", default=4, type=int, 
                    	help="specify the number of time spots, default: 4,                     \n \
                          must be provided and valid. ")
  parser.add_argument("-t", "--transFunc", dest="transFunc", default='simple', choices=['simple', 'SD'],
                      help="specify the method to summarize replicates data, default: simple, \n \
                          choices: simple, SD                                                 \n \
                          NOTE:                                                               \n \
                          simple: simple averaging                                            \n \
                          SD: standard deviation weighted averaging;" )
  parser.add_argument("-f", "--fillMethod", dest="fillMethod", default='none', choices=['none', 'zero', 'linear', 'quadratic', 'cubic', 'slinear', 'nearest'],
                    	help= "specify the method to fill missing, default: none,               \n \
                          choices: none, zero, linear, quadratic, cubic, slinear, nearest  \n \
                          NOTE:                                            \n \
                          none: fill up with zeros;                        \n \
                          zero: fill up with zero order splines;           \n \
                          linear: fill up with linear splines;             \n \
                          slinear: fill up with slinear;                   \n \
                          quadratic: fill up with quadratic spline;             \n \
                          cubic: fill up with cubic spline;                \n \
                          nearest: fill up with nearest neighbor") 
  
  arg_namespace = parser.parse_args()

  #get arguments
  print >>sys.stderr, "lsa-compute"
  print >>sys.stderr, "copyright Li Xia, lxia@usc.edu"
  print >>sys.stderr, "learning arguments..."
  
  delayLimit = vars(arg_namespace)['delayLimit']
  fillMethod = vars(arg_namespace)['fillMethod']
  permuNum = vars(arg_namespace)['permuNum']
  dataFile = vars(arg_namespace)['dataFile']				#dataFile
  resultFile = vars(arg_namespace)['resultFile']			#resultFile
  repNum = vars(arg_namespace)['repNum']
  transFunc = vars(arg_namespace)['transFunc']
  bootNum = vars(arg_namespace)['bootNum']
  spotNum = vars(arg_namespace)['spotNum']

  #check transFunc and repNum compatibility
  if repNum <= 5 and transFunc == 'SD':
    print >>sys.stderr, "Not enough replicates for SD-weighted averaging, fall back to simple"
    transFunc = 'simple'
  if transFunc == 'simple':
    fTransform = lsalib.simpleAverage
  else:
    fTransform = lsalib.sdAverage
  #print "\t".join([ 'delayLimit', 'fillMethod', 'permuNum', 'dataFile', 'resultFile', 'repNum', 'spotNum', 'bootNum', 'transFunc' ])
  #print "\t".join(['%s']*9) % (delayLimit, fillMethod, permuNum, dataFile, resultFile, repNum, spotNum, bootNum, transFunc)
  
  #start timing main
  start_time = time.time()

  #input
  try:
    print >>sys.stderr, "reading raw data from dataFile..."
    rawData=np.genfromtxt( dataFile, comments='#', delimiter='\t', missing_values=['na',''], filling_values=np.nan, usecols=range(1,spotNum*repNum+1) )
    dataFile.seek(0)  #rewind
    factorLabels=list(np.genfromtxt( dataFile, comments='#', delimiter='\t', usecols=xrange(0,1), dtype='string' ))
  except:
    print >>sys.stderr, "unexpected error:", sys.exc_info()[0]
    print >>sys.stderr, "error reading dataFile, please check the input format, spotNum and repNum \n \
                         input shall be a tab delimited txt file with '#' lines as comments and first column as factor label. \n \
                         After that, it shall have spotNum * repNum numeric cells for repNum-replicated spotNum-spotted series data. "
    exit(0)

  ###print rawData, factorLabels
  factorNum = rawData.shape[0]
  cleanData=np.zeros( ( factorNum, repNum, spotNum), dtype='float' ) # (num_rows-1) x (num_cols/repNum) x (repNum)
  for i in xrange(0, factorNum):
    for j in xrange(0, repNum):
      #print rawData[i], j, spotNum*repNum, repNum, np.arange(j,spotNum*repNum, repNum)
      cleanData[i,j] = rawData[i][np.arange(j,spotNum*repNum,repNum)]
  ###print cleanData

  for i in xrange(0, factorNum):
    for j in range(0, repNum):
      cleanData[i,j] = lsalib.fillMissing( cleanData[i,j], fillMethod )
    
  #calculation
  #[ Seq X's Idx, Seq Y's Idx, LS Score, CI_low, CI_high, X's Start Position, 
  #        Y's Start Position, Alignment Length, X delay to Y,
  #        P-value, Pearson' Correlation, P-value of PCC, Q-value ]
  print >>sys.stderr, "data size factorNum, repNum, spotNum = %s, %s, %s" % (cleanData.shape[0], cleanData.shape[1], cleanData.shape[2])
  print >>sys.stderr, "calculating ..."
  lsaTable = lsalib.applyAnalysis( cleanData, delayLimit=delayLimit, bootNum=bootNum, permuNum=permuNum, fTransform=fTransform )
  print >>sys.stderr, "writing results ..."
  print >>resultFile,  "\t".join(['X','Y','LS','lowCI','upCI','Xs','Ys','Len','Delay','P','PCC','Ppcc','Q'])
  for row in lsaTable:
    #print [factorLabels[row[0]], factorLabels[row[1]]] + ["%.4f" % v if isinstance(v, float) else v for v in row[2:]]
    print >>resultFile, "\t".join(['%s']*13) % \
          tuple([factorLabels[row[0]], factorLabels[row[1]] ] + ["%.4f" % v if isinstance(v, float) else v for v in row[2:]])

  #output

  print >>sys.stderr, "finishing up..."
  resultFile.close()
  end_time=time.time()
  print >>sys.stderr, "time elapsed %f seconds" % (end_time-start_time)
  print >>sys.stderr, "Thank you for using lsa-compute, byebye"

if __name__=="__main__":
  main()
