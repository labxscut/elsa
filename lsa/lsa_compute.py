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
  #debug import
  from . import lsalib
except ImportError:
  #install import
  from lsa import lsalib
  #np.seterr(all='raise')
import lsa

# let's display something about VERSION first, need to put VERSION.py there
#execfile(os.path.join(os.path.dirname(os.path.dirname(lsa.__file__)),\
#			    'VERSION.py')) #only one variable named version_desc

def main():  

  __script__ = "lsa_compute"
  version_desc = lsalib.safeCmd('lsa_version')
  version_print = "%s (rev: %s) - copyright Li Charlie Xia, lcxia@scut.edu.cn" \
    % (__script__, version_desc) 
  print(version_print, file=sys.stderr)

  #version_desc = "lsa_compute (rev: %s) - copyright Li Charlie Xia, lixia@stanford.edu" \
  #    % open(os.path.join(os.path.dirname(os.path.dirname(lsa.__file__)), 'VERSION.txt')).read().strip()
  #lsa.__version__
  #print >>sys.stderr, version_desc
  #print >>sys.stderr, "lsa package revision and tag:", 
  #print >>sys.stderr, "for package version, see VERSION.txt in your lsa package installation path"

  # define arguments: delayLimit, fillMethod, pvalueMethod
  parser = argparse.ArgumentParser()

  arg_precision_default=1000
  arg_delayLimit_default=0

  parser.add_argument("dataFile", metavar="dataFile", type=argparse.FileType('r'), \
      help="the input data file,\n \
      m by (r * s)tab delimited text; top left cell start with \
      '#' to mark this is the header line; \n \
      m is number of variables, r is number of replicates, \
      s it number of time spots; \n \
      first row: #header  s1r1 s1r2 s2r1 s2r2; \
      second row: x  ?.?? ?.?? ?.?? ?.??; for a 1 by (2*2) data")
  parser.add_argument("resultFile", metavar="resultFile", type=argparse.FileType('w'), \
      help="the output result file")
  parser.add_argument("-e", "--extraFile", dest="extraFile", default=None, \
      type=argparse.FileType('r'),
      help="specify an extra datafile, otherwise the first datafile will be used \n \
            and only lower triangle entries of pairwise matrix will be computed")
  parser.add_argument("-d", "--delayLimit", dest="delayLimit", default=arg_delayLimit_default, type=int,\
      help="specify the maximum delay possible, default: {},\n \
            must be an integer >=0 and <spotNum".format(arg_delayLimit_default))
  parser.add_argument("-m", "--minOccur", dest="minOccur", default=50, type=int, 
      help="specify the minimum occurence percentile of all times, default: 50,\n")
  parser.add_argument("-p", "--pvalueMethod", dest="pvalueMethod", default="perm", \
      choices=["perm", "theo", "mix"],
      help="specify the method for p-value estimation, \n \
            default: pvalueMethod=perm, i.e. use  permutation \n \
            theo: theoretical approximaton; if used also set -a value. \n \
            mix: use theoretical approximation for pre-screening \
            if promising (<0.05) then use permutation. ")
  parser.add_argument("-x", "--precision", dest="precision", default=arg_precision_default, type=int,\
      help="permutation/precision, specify the permutation \n \
            number or precision=1/permutation for p-value estimation. \n \
            default is {}, must be an integer >0 ".format(arg_precision_default) )
  parser.add_argument("-b", "--bootNum", dest="bootNum", default=0, type=int, \
      choices=[0, 100, 200, 500, 1000, 2000],
      help="specify the number of bootstraps for 95%% confidence \
            interval estimation, default: 100,\n \
            choices: 0, 100, 200, 500, 1000, 2000. \n \
            Setting bootNum=0 avoids bootstrap. \n \
            Bootstrap is not suitable for non-replicated data.")   #use %% to print %
  parser.add_argument("-r", "--repNum", dest="repNum", default=1, type=int,
      help="specify the number of replicates each time spot, default: 1,\n \
            must be provided and valid. ")
  parser.add_argument("-s", "--spotNum", dest="spotNum", default=4, type=int, 
      help="specify the number of time spots, default: 4,\n \
            must be provided and valid. ")
  parser.add_argument("-t", "--transFunc", dest="transFunc", default='simple', \
      choices=['simple', 'SD', 'Med', 'MAD'],\
      help="specify the method to summarize replicates data, default: simple, \n \
            choices: simple, SD, Med, MAD                                     \n \
            NOTE:                                                             \n \
            simple: simple averaging                                          \n \
            SD: standard deviation weighted averaging                         \n \
            Med: simple Median                                                \n \
            MAD: median absolute deviation weighted median;" )
  parser.add_argument("-f", "--fillMethod", dest="fillMethod", default='none', \
      choices=['none', 'zero', 'linear', 'quadratic', 'cubic', 'slinear', 'nearest'], \
      help="specify the method to fill missing, default: none,               \n \
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
  parser.add_argument("-n", "--normMethod", dest="normMethod", default='robustZ', \
      choices=['percentile', 'percentileZ', 'pnz', 'robustZ', 'rnz', 'none'], \
      help="must specify the method to normalize data, default: robustZ, \n \
            choices: percentile, none, pnz, percentileZ, robustZ or a float  \n \
            NOTE:                                                   \n \
            percentile: percentile normalization, including zeros (only with perm)\n \
            pnz: percentile normalization, excluding zeros (only with perm) \n  \
            percentileZ: percentile normalization + Z-normalization \n \
            rnz: percentileZ normalization + excluding zeros + robust estimates (theo, mix, perm OK) \n \
            robustZ: percentileZ normalization + robust estimates \n \
            (with perm, mix and theo, and must use this for theo and mix, default) \n")
  parser.add_argument("-q", "--qvalueMethod", dest="qvalueMethod", \
      default='scipy', choices=['scipy'],
      help="specify the qvalue calculation method, \n \
            scipy: use scipy and storeyQvalue function, default \n \
            ")
            #R: use R's qvalue package, require X connection")
  parser.add_argument("-T", "--trendThresh", dest="trendThresh", default=None, \
      type=float, \
      help="if trend series based analysis is desired, use this option \n \
            NOTE: when this is used, must also supply reasonble \n \
            values for -p, -a, -n options")
  parser.add_argument("-a", "--approxVar", dest="approxVar", default=1, type=float,\
      help="if use -p theo and -T, must set this value appropriately, \n \
            precalculated -a {1.25, 0.93, 0.56,0.13 } for i.i.d. standard normal null \n \
            and -T {0, 0.5, 1, 2} respectively. For other distribution \n \
            and -T values, see FAQ and Xia et al. 2013 in reference")
  parser.add_argument("-v", "--progressive", dest="progressive", default=0, type=int, 
      help="specify the number of progressive output to save memory, default: 0,\n \
            2G memory is required for 1M pairwise comparison. ")
  arg_namespace = parser.parse_args()
  
  delayLimit = vars(arg_namespace)['delayLimit']
  fillMethod = vars(arg_namespace)['fillMethod']
  normMethod = vars(arg_namespace)['normMethod']
  qvalueMethod = vars(arg_namespace)['qvalueMethod']
  pvalueMethod = vars(arg_namespace)['pvalueMethod']
  precision = vars(arg_namespace)['precision']
  minOccur = vars(arg_namespace)['minOccur']
  dataFile = vars(arg_namespace)['dataFile']				#dataFile
  extraFile = vars(arg_namespace)['extraFile']				#extraFile
  resultFile = vars(arg_namespace)['resultFile']			#resultFile
  repNum = vars(arg_namespace)['repNum']
  transFunc = vars(arg_namespace)['transFunc']
  bootNum = vars(arg_namespace)['bootNum']
  spotNum = vars(arg_namespace)['spotNum']
  approxVar = vars(arg_namespace)['approxVar'] 
  trendThresh = vars(arg_namespace)['trendThresh'] 
  progressive = vars(arg_namespace)['progressive'] 

  try:
    extraFile_name = extraFile.name 
  except AttributeError:
    extraFile_name = ''

  #determin if trend series analysis and see if trendThresh option is valid
  assert trendThresh==None or trendThresh>=0

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
    print("Not enough replicates for SD-weighted averaging, fall back to simpleAverage", file=sys.stderr)
    transFunc = 'simple'

  if repNum < 5 and ( transFunc == 'MAD' ):
    print("Not enough replicates for Median Absolute Deviation, fall back to simpleMedian", file=sys.stderr)
    transFunc = 'Med'

  #check normMethod
  if normMethod == 'none':
    zNormalize = lsalib.noneNormalize
  elif normMethod == 'percentile':
    zNormalize = lsalib.percentileNormalize
  elif normMethod == 'percentileZ':
    zNormalize = lsalib.percentileZNormalize
  elif normMethod == 'robustZ':
    zNormalize = lsalib.robustZNormalize
  elif normMethod == 'pnz':
    zNormalize = lsalib.noZeroNormalize
  elif normMethod == 'rnz':
    zNormalize = lsalib.robustNoZeroNormalize
  else:
    zNormalize = lsalib.percentileZNormalize # fallback to default

  assert precision>0, "precision %s is not positive" % str(precision) 
  
  print("\t".join(['delayLimit','minOccur','fillMethod','pvalueMethod',\
      'precision','dataFile','extraFile','resultFile','repNum','spotNum',\
      'bootNum','transFunc','normMethod','approxVar','trendThresh']))
  print("\t".join(['%s']*15) % (delayLimit,minOccur,fillMethod,pvalueMethod,\
      precision,dataFile.name,extraFile_name,resultFile.name,repNum,spotNum,\
      bootNum,transFunc,normMethod,str(approxVar),str(trendThresh)))
  
  #start timing main
  start_time = time.time()
  #print("starting1...")

  #datafile handling
  onDiag = False
  try:
    firstData=np.genfromtxt( \
        dataFile, comments='#', delimiter='\t', missing_values=['na','','NA'], \
        filling_values=np.nan, usecols=list(range(1,spotNum*repNum+1)) )
    #print("starting10...")
    if len(firstData.shape)==1:
      firstData=np.array([firstData])
    #print("starting101...")
    #print firstData.shape
    dataFile.seek(0)  #rewind
    #print("starting102...")
    #print "we can get here"
    firstFactorLabels=np.genfromtxt( \
        dataFile, comments='#', delimiter='\t', \
        usecols=range(0,1), dtype='str' )
    #print(firstFactorLabels)
    firstFactorLabels=firstFactorLabels.tolist()
    #print("starting103...")
    if type(firstFactorLabels)==str:
      firstFactorLabels=[firstFactorLabels]
    #print("starting11...")
    #print "but we cann't get here"
    if not extraFile:
      #print("starting12...")
      onDiag = True
      #print >>sys.stderr, "reading raw data from dataFile..."
      dataFile.seek(0)  #rewind
      secondData=np.genfromtxt( dataFile, comments='#', delimiter='\t', \
          missing_values=['na','','NA'], \
          filling_values=np.nan, usecols=list(range(1,spotNum*repNum+1)) )
      if len(secondData.shape)==1:
        secondData=np.array([secondData.shape])
      dataFile.seek(0)  #rewind
      secondFactorLabels=np.genfromtxt( dataFile, comments='#', delimiter='\t', \
          usecols=range(0,1), dtype='str' ).tolist()
      if type(secondFactorLabels)==str:
        secondFactorLabels=[secondFactorLabels]
    else:
      #print("starting13...")
      #print "can we get here"
      #dataFile.close()  #incase feeding the same file twice
      extraFile.seek(0)
      secondData=np.genfromtxt( \
          extraFile, comments='#', delimiter='\t', missing_values=['na','','NA'], \
          filling_values=np.nan, usecols=list(range(1,spotNum*repNum+1)) )
      if len(secondData.shape)==1:
        secondData=np.array([secondData.shape])
      print(secondData.shape)
      extraFile.seek(0)  #rewind
      secondFactorLabels=np.genfromtxt( extraFile, comments='#', delimiter='\t', \
          usecols=range(0,1), dtype='str' ).tolist()
      if type(secondFactorLabels)==str:
        secondFactorLabels=[secondFactorLabels]
  except ValueError:
    print("ValueError:", str(ValueError), file=sys.stderr)
    exit(0)
  except TypeError:
    print("TypeError:", str(TypeError), file=sys.stderr)
    exit(0)
  except:
    print("unexpected error:", sys.exc_info()[0], file=sys.stderr)
    print("error reading dataFile, \n \
    please check the input format, spotNum and repNum \n \
    input shall be a tab delimited txt file with first \n \
    line starts with '#' as column names and first column as factor labeles. \n \
    it allows other rows start with '#' to be comment lines \n \
    An in total, it shall have spotNum * repNum numeric cells \n \
    for repNum-replicated spotNum-timepoint series data. ", file=sys.stderr)
    exit(0)
  #print("starting2...")

  ###print rawData, factorLabels
  cleanData = []
  for rawData in [firstData, secondData]:
    factorNum = rawData.shape[0]
    tempData=np.zeros( ( factorNum, repNum, spotNum), dtype='float' ) 
    # (num_rows-1) x (num_cols/repNum) x (repNum)
    for i in range(0, factorNum):
      for j in range(0, repNum):
        try:
          tempData[i,j] = rawData[i][np.arange(j,spotNum*repNum,repNum)]
        except IndexError:
          print("""Error: one input file need more than two data row
or use -e to specify another input file""", file=sys.stderr)
          quit()
    for i in range(0, factorNum):
      for j in range(0, repNum):
        tempData[i,j] = lsalib.fillMissing( tempData[i,j], fillMethod )
    cleanData.append(tempData)
  #print tempData
  #print("starting3...")
    
  #calculation
  #[ Seq X's Idx, Seq Y's Idx, LS Score, CI_low, CI_high, X's Start Position, 
  #        Y's Start Position, Alignment Length, X delay to Y,
  #        P-value, Pearson' Correlation, P-value of PCC, Q-value ]
  print("firstData factorNum, repNum, spotNum = %s, %s, %s" \
      % (cleanData[0].shape[0], cleanData[0].shape[1], cleanData[0].shape[2]), file=sys.stderr)
  print("secondData factorNum, repNum, spotNum = %s, %s, %s" \
      % (cleanData[1].shape[0], cleanData[1].shape[1], cleanData[1].shape[2]), file=sys.stderr)
  print("calculating ", cleanData[0].shape[0]*cleanData[1].shape[0], "pairwise local similarity scores...", file=sys.stderr)
  lsalib.applyAnalysis(cleanData[0], cleanData[1], onDiag=onDiag, \
      delayLimit=delayLimit, bootNum=bootNum, minOccur=minOccur/100.,\
      pvalueMethod=pvalueMethod, precisionP=precision, fTransform=fTransform,\
      zNormalize=zNormalize, approxVar=approxVar, resultFile=resultFile,\
      firstFactorLabels=firstFactorLabels, trendThresh=trendThresh,\
      secondFactorLabels=secondFactorLabels, qvalueMethod=qvalueMethod, progressive=progressive)

  #print >>sys.stderr, "writing results ..."
  #col_labels= ['X','Y','LS','lowCI','upCI','Xs','Ys','Len','Delay','P','PCC','Ppcc','SPCC','Pspcc','SCC','Pscc','SSCC','Psscc',
  #    'Q','Qpcc','Qspcc','Qscc','Qsscc']
  #print >>resultFile,  "\t".join(col_labels)

  #print lsaTable
  #for row in lsaTable:
    #print [factorLabels[row[0]], factorLabels[row[1]]] + ["%.4f" % v if isinstance(v, float) else v for v in row[2:13]]
  #  print >>resultFile, "\t".join(['%s']*len(col_labels)) % \
  #    tuple([firstFactorLabels[row[0]], secondFactorLabels[row[1]] ] + ["%.4f" % np.round(v, decimals=4) if isinstance(v, float) else v for v in row[2:]])

  print("finishing up...", file=sys.stderr)
  resultFile.close()
  end_time=time.time()
  print("time elapsed %f seconds" % (end_time-start_time), file=sys.stderr)
  #print >>sys.stderr, "Thank you for using lsa-compute, byebye"

if __name__=="__main__":
  #print >>sys.stderr, lsaver.rev
  main()
