#!/usr/bin/env python
#la_compute -- computation script for Liquid Association calculation 

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

def main():  

  # define arguments: delayLimit, fillMethod, pvalueMethod
  parser = argparse.ArgumentParser(description="New Liquid Association Commandline Tool")

  parser.add_argument("dataFile", metavar="dataFile", type=argparse.FileType('r'), help="the input data file,\n \
                        m by (r * s)tab delimited text; top left cell start with '#' to mark this is the header line; \n \
                        m is number of variables, r is number of replicates, s it number of time spots; \n \
                        first row: #header  s1r1 s1r2 s2r1 s2r2; second row: x  ?.?? ?.?? ?.?? ?.??; for a 1 by (2*2) data")
  parser.add_argument("scoutFile", metavar="scoutFile", type=argparse.FileType('r'),
                        help="the input datafile specify the scouting pairs, \n \
                            it can be any tab delimited file (e.g. .lsa) with (xi, yi) pair indecies for scouting pairs")
  parser.add_argument("resultFile", metavar="resultFile", type=argparse.FileType('w'), help="the output result file")
  parser.add_argument("-x", "--xiCol", dest="xiCol", default=26, type=int, 
                    	help="specify the x-th column to store Xi indecies")
  parser.add_argument("-y", "--yiCol", dest="yiCol", default=27, type=int, 
                    	help="specify the y-th column to store Yi indecies")
  parser.add_argument("-p", "--pvalueMethod", dest="pvalueMethod", default=1000, type=int,
                    	help="specify the method=sgn(pvalueMethod) and precision=1/abs(pvalueMethod) for p-value estimation, \n \
                          default: pvalueMethod=1000, i.e. precision=0.001 and mode=permutation \n \
                          mode +: permutation approximaton; -: unused, fall back to permutation")
  parser.add_argument("-m", "--minOccur", dest="minOccur", default=50, type=int,
                              help="specify the minimum occurence percentile of all times, default: 50,\n")
  parser.add_argument("-b", "--bootNum", dest="bootNum", default=0, type=int, choices=[0, 100, 200, 500, 1000, 2000],
                    	help="specify the number of bootstraps for 95%% confidence interval estimation, default: 100,\n \
                          choices: 0, 100, 200, 500, 1000, 2000. \n \
                          Setting bootNum=0 avoids bootstrap. Bootstrap is not suitable for non-replicated data.")   #use %% to print %
  parser.add_argument("-r", "--repNum", dest="repNum", default=1, type=int,
                    	help="specify the number of replicates each time spot, default: 1,      \n \
                          must be provided and valid. ")
  parser.add_argument("-s", "--spotNum", dest="spotNum", default=4, type=int, 
                    	help="specify the number of time spots, default: 4,                     \n \
                          must be provided and valid. ")
  parser.add_argument("-t", "--transFunc", dest="transFunc", default='simple', choices=['simple', 'SD', 'Med', 'MAD'],
                        help="specify the method to summarize replicates data, default: simple, \n \
                          choices: simple, SD, Med, MAD                                       \n \
                          NOTE:                                                               \n \
                          simple: simple averaging                                            \n \
                          SD: standard deviation weighted averaging                           \n \
                          Med: simple Median                                                  \n \
                          MAD: median absolute deviation weighted median;" )
  parser.add_argument("-f", "--fillMethod", dest="fillMethod", default='linear', 
                        choices=['none', 'zero', 'linear', 'quadratic', 'cubic', 'slinear', 'nearest'],
                    	help= "specify the method to fill missing, default: none,               \n \
                          choices: none, zero, linear, quadratic, cubic, slinear, nearest  \n \
                          operation AFTER normalization:  \n \
                          none: fill up with zeros ;   \n \
                          operation BEFORE normalization:  \n \
                          zero: fill up with zero order splines;           \n \
                          linear: fill up with linear splines;             \n \
                          slinear: fill up with slinear;                   \n \
                          quadratic: fill up with quadratic spline;             \n \
                          cubic: fill up with cubic spline;                \n \
                          nearest: fill up with nearest neighbor") 
  parser.add_argument("-n", "--normMethod", dest="normMethod", default='pnz', #choices=['percentile', 'percentileZ', 'pnz', 'none'],
                        help= "specify the method to normalize data, default: percentile,       \n \
                          choices: percentile, none, pnz, percentileZ                         \n \
                          NOTE:                                                         \n \
                          percentile: percentile normalization, including zeros                                    \n \
                          pnz: percentile normalization, excluding zeros                                    \n \
                          percentileZ: percentile normalization + Z-normalization                                   \n \
                          none or a float number for variance: no normalization and calculate Ptheo with user specified variance, default=1 ")
  
  arg_namespace = parser.parse_args()

  #get arguments
  print >>sys.stderr, "la_compute ($Revision$) - copyright Li Charlie Xia, lxia@usc.edu"
  
  dataFile = vars(arg_namespace)['dataFile']				#dataFile
  scoutFile = vars(arg_namespace)['scoutFile']				#extraFile
  resultFile = vars(arg_namespace)['resultFile']			#resultFile
  xiCol = vars(arg_namespace)['xiCol']
  yiCol = vars(arg_namespace)['yiCol']
  minOccur = vars(arg_namespace)['minOccur']
  pvalueMethod = vars(arg_namespace)['pvalueMethod']
  bootNum = vars(arg_namespace)['bootNum']
  repNum = vars(arg_namespace)['repNum']
  spotNum = vars(arg_namespace)['spotNum']
  transFunc = vars(arg_namespace)['transFunc']
  fillMethod = vars(arg_namespace)['fillMethod']
  normMethod = vars(arg_namespace)['normMethod']

  #assign transform function
  if transFunc == 'SD':
    fTransform = lsalib.sdAverage
  elif transFunc == 'Med':
    fTransform = lsalib.simpleMedian   # Median
  elif transFunc == 'MAD':
    fTransform = lsalib.madMedian      # Median/MAD
  else:
    fTransform = lsalib.simpleAverage   # fallback to default Avg
  
  #check transFunc and repNum compatibility
  if repNum < 5 and ( transFunc == 'SD' ):
    print >>sys.stderr, "Not enough replicates for SD-weighted averaging, fall back to simpleAverage"
    transFunc = 'simple'

  if repNum < 5 and ( transFunc == 'MAD' ):
    print >>sys.stderr, "Not enough replicates for Median Absolute Deviation, fall back to simpleMedian"
    transFunc = 'Med'

  #check normMethod
  if normMethod == 'none':
    zNormalize = lsalib.noneNormalize
  elif normMethod == 'percentile':
    zNormalize = lsalib.percentileNormalize
  elif normMethod == 'pnz':
    zNormalize = lsalib.noZeroNormalize  
  else:
    zNormalize = lsalib.noZeroNormalize # fallback to default
  
  pars = ['fillMethod','minOccur', 'pvalueMethod','dataFile','scoutFile','resultFile','repNum','spotNum','bootNum','transFunc','normMethod']
  print "\t".join(pars)
  print "\t".join(['%s']*len(pars)) \
      % (fillMethod,minOccur,pvalueMethod,dataFile.name,scoutFile.name,resultFile.name,repNum,spotNum,bootNum,transFunc,normMethod)
  
  #start timing main
  start_time = time.time()

  #datafile handling
  try:
    inputData=np.genfromtxt( dataFile, comments='#', delimiter='\t', missing_values=['na','','NA'], \
        filling_values=np.nan, usecols=range(1,spotNum*repNum+1) )
    dataFile.seek(0)  #rewind
    factorLabels=np.genfromtxt( dataFile, comments='#', delimiter='\t', usecols=xrange(0,1), dtype='string' )
    if factorLabels.shape != ():
      factorLabels=list(factorLabels)
    else: #in case of single row
      factorLabels=[str(factorLabels)]
      inputData.shape=(1,inputData.shape[0])
  except:
    dataFile.seek(0)
    print >>sys.stderr, [str(np.genfromtxt( dataFile, comments='#', delimiter='\t', usecols=xrange(0,1), dtype='string' ))]
    print >>sys.stderr, spotNum*repNum+1
    print >>sys.stderr, "unexpected error:", sys.exc_info()[0]
    print >>sys.stderr, "error reading dataFile, please check the input format, spotNum and repNum \n \
                         input shall be a tab delimited txt file with first line starts with '#' as column names and first column as factor labeles. \n \
                         it allows other rows start with '#' to be comment lines \n \
                         An in total, it shall have spotNum * repNum numeric cells for repNum-replicated spotNum-timepoint series data. "
    exit(0)

  #print inputData.shape
  factorNum = inputData.shape[0]
  tempData=np.zeros( (factorNum, repNum, spotNum), dtype='float' ) # (num_rows-1) x (num_cols/repNum) x (repNum)
  for i in xrange(0, factorNum):
    for j in xrange(0, repNum):
      tempData[i,j] = inputData[i][np.arange(j,spotNum*repNum,repNum)]
  for i in xrange(0, factorNum):
    for j in range(0, repNum):
      tempData[i,j] = lsalib.fillMissing( tempData[i,j], fillMethod )
    
  scoutVars = []
  ri = 0
  for row in csv.reader( scoutFile, delimiter="\t" ):
    if ri == 0:
      ri = ri + 1 #skip first row
      continue
    scoutVars.append( [int(row[xiCol-1]),int(row[yiCol-1])] )

  #print tempData
  #print scoutVars
  #print factorLabels

  lsalib.applyLA(tempData, scoutVars, factorLabels, bootNum=bootNum, minOccur=minOccur/100.,\
      pvalueMethod=pvalueMethod, fTransform=fTransform, zNormalize=zNormalize, resultFile=resultFile)

  print >>sys.stderr, "finishing up..."
  resultFile.close()
  end_time=time.time()
  print >>sys.stderr, "time elapsed %f seconds" % (end_time-start_time)

if __name__=="__main__":
  main()
