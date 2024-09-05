#lalib.py -- Library of Liquid Association Analysis(LAA) Package
#LICENSE: BSD

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

"""lalib.py -- Library of Liquid Association Analysis(LAA) Package

  NOTE: numpy and scipy is required to use thie module
  NOTE: accepts input sequence table as delimited text file 
        with first row is the "Date" and other factor labels and
        first column is the time spot labels
"""

import csv, sys, os, random
import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.stats
#R through Rpy
rpy_import=False
#try:
#  import rpy2.rlike.container as rlc
#  import rpy2.robjects as ro
#  from rpy2.robjects.numpy2ri import numpy2ri
#  ro.conversion.py2ri = numpy2ri
#  r = ro.r
#
#  #print '''setwd("%s")''' % os.environ.get('PWD')
#  r('''setwd("%s")''' % os.environ.get('PWD'))
#  r('''options(warn=-1)''')
#except ImportError:
#  rpy_import=False
#  print >>sys.stderr, "IMPORTANT!!!: R and rpy2 are not working on this system"
#  print >>sys.stderr, "IMPORTANT!!!: All calculatios fall back to scipy"
#
#import lower level resource
try:
  #else run as installed
  from lsa import compcore
  from lsa import lsalib
  np.seterr(all='ignore')                             #ignore RuntimeWarning of abs in installed mode
except ImportError:
  #try for debug
  from . import compcore
  from . import lsalib
  np.seterr(all='warn')

#global variable, stores calculated p-values.
#P_table = dict()
disp_decimal=8
kcut_min=100
Rmax_min=10
my_decimal = 2    # preset x step size for P_table
pipi = np.pi**2 # pi^2
pipi_inv = 1/pipi
Q_lam_step = 0.05
Q_lam_max = 0.95

### LAA functions ###

def applyLA(inputData, scoutVars, factorLabels, bootCI=.95, bootNum=1000, minOccur=.50, \
    pvalueMethod="perm", precision=1000,\
    fTransform=lsalib.simpleAverage, zNormalize=lsalib.noZeroNormalize, resultFile=None, qvalue_func=lsalib.storeyQvalue):

  col_labels = ['X','Y','Z','LA','lowCI','upCI','P','Q','Xi','Yi','Zi']
  print("\t".join(col_labels), file=resultFile)

  #print(inputData.shape)
  inputFactorNum = inputData.shape[0]
  inputRepNum = inputData.shape[1]
  inputSpotNum = inputData.shape[2]
  scoutNum = len(scoutVars)
  cp = np.array( [False]*inputFactorNum**3, dtype='bool' ) #consider bitvector
  cp.shape = (inputFactorNum, inputFactorNum, inputFactorNum)
  laTable = [None]*inputFactorNum*scoutNum
  pvalues = np.zeros(inputFactorNum*scoutNum, dtype='float')
  timespots = inputSpotNum #same length already assumed
  replicates = inputRepNum
  ti = 0
  for i in range(0, scoutNum):
    Xi = scoutVars[i][0] - 1
    Yi = scoutVars[i][1] - 1
    #print Xi, Yi
    Xo = np.ma.masked_invalid(inputData[Xi], copy=True) #need to convert to masked array with na's, not F-normalized
    Yo = np.ma.masked_invalid(inputData[Yi], copy=True) #need to convert to masked array with na's, not F-normalized
    for j in range(0, inputFactorNum):
      Zi = j
      if Xi == Yi or Xi == Zi or Zi == Yi:
        continue   #ignore invalid entries
      if cp[Xi,Yi,Zi] or cp[Xi,Zi,Yi] or cp[Zi,Xi,Yi] or cp[Zi,Yi,Xi] or cp[Yi,Xi,Zi] or cp[Yi,Zi,Xi]:
        continue   #ignore redundant entries
      cp[Xi,Yi,Zi]=True
      Zo = np.ma.masked_invalid(inputData[Zi], copy=True)    # need to convert to masked array with na's, not F-normalized
      Xo_badOccur = np.sum(np.logical_not(np.isnan(lsalib.ma_average(Xo)), lsalib.ma_average(Xo)==0))/float(timespots) < minOccur
      Yo_badOccur = np.sum(np.logical_not(np.isnan(lsalib.ma_average(Yo)), lsalib.ma_average(Yo)==0))/float(timespots) < minOccur
      Zo_badOccur = np.sum(np.logical_not(np.isnan(lsalib.ma_average(Zo)), lsalib.ma_average(Zo)==0))/float(timespots) < minOccur
      if Xo_badOccur or Yo_badOccur or Zo_badOccur:           #either one of these not satisfying the minOccur criteria
        continue						 # minOccur not satisfied for one of the factors
      if np.all(Xo.mask) or np.all(Yo.mask) or np.all(Zo.mask):  # not any unmasked value in Xz or Yz, all nan in input, continue
        continue
      LA_score = singleLA(Xo, Yo, Zo, fTransform, zNormalize)                          # do LA computation

      if np.isnan(LA_score):
        print("found na in LA score, to be fixed", file=sys.stderr)
        #print >>sys.stderr, Xo, zNormalize(fTransform(Xo))
        #print >>sys.stderr, Yo, zNormalize(fTransform(Yo))
        #print >>sys.stderr, Zo, zNormalize(fTransform(Zo))
        laTable[ti] = [Xi, Yi, Zi, LA_score, LA_score, LA_score, LA_score]
        ti += 1
        continue

      #This Part Must Follow Static Calculation Part
      #Otherwise Yz may be changed, now it is copied
      #np.ma.array(copy=True) to copy, otherwise is only reference
      if pvalueMethod == "perm":
        Xp = np.ma.array(Xo,copy=True)
        Yp = np.ma.array(Yo,copy=True)
        Zp = np.ma.array(Zo,copy=True)
        laP = LApermuPvalue(Xp, Yp, Zp, precision, np.abs(LA_score), fTransform, zNormalize)          # do Permutation Test
      else: #reserved 
        print("this branch should not be receahed by now", file=sys.stderr)

      pvalues[ti] = laP
      #print "bootstrap computing..."
      if bootNum > 0: #do BS
        Xb = np.ma.array(Xo,copy=True)
        Yb = np.ma.array(Yo,copy=True)
        Zb = np.ma.array(Zo,copy=True)
        (LA_score, Sl, Su) = LAbootstrapCI(Xb, Yb, Zb, LA_score, bootCI, bootNum, fTransform, zNormalize)        # do Bootstrap CI
      else: #skip BS
        (LA_score, Sl, Su) = (LA_score, LA_score, LA_score)

      laTable[ti] = [Xi, Yi, Zi, LA_score, Sl, Su, laP]
      ti += 1

  pvalues = pvalues[:ti]
  laTable = laTable[:ti]
  qvalues = qvalue_func( pvalues )
  for k in range(0, len(qvalues)):
    laTable[k] = laTable[k] + [ qvalues[k], laTable[k][0]+1, laTable[k][1]+1, laTable[k][2]+1 ]

  #print laTable
  for row in laTable:
    print("\t".join( ['%s']*len(col_labels) ) % \
      tuple( [factorLabels[row[0]], factorLabels[row[1]], factorLabels[row[2]]] \
          + ["%f" % np.round(v, decimals=disp_decimal) if isinstance(v, float) else v for v in row[3:]] ), file=resultFile)

def singleLA(series1, series2, series3, fTransform, zNormalize):
  return compcore.calc_LA(zNormalize(fTransform(series1)),zNormalize(fTransform(series2)),zNormalize(fTransform(series3)))

def calc_LA(series1, series2, series3):
  n1 = len(series1)
  n2 = len(series2)
  n3 = len(series3)
  assert n1==n2 and n2 == n3
  return np.sum(series1*series2*series3)/n1

def LAbootstrapCI(series1, series2, series3, LA_score, bootCI, bootNum, fTransform, zNormalize, debug=0):
  ### no feasible, skipping bootstraping
  if series1.shape[0] == 1:
    return (LA_score, LA_score, LA_score)

  BS_set = np.zeros(bootNum, dtype='float')
  for i in range(0, bootNum):
    Xb = np.ma.array([ sample_wr(series1[:,j], series1.shape[0]) for j in range(0,series1.shape[1]) ]).T
    Yb = np.ma.array([ sample_wr(series2[:,j], series2.shape[0]) for j in range(0,series2.shape[1]) ]).T
    Zb = np.ma.array([ sample_wr(series3[:,j], series3.shape[0]) for j in range(0,series3.shape[1]) ]).T
    BS_set[i] = compcore.calc_LA(Xb, Yb, Zb)
  BS_set.sort()                                 #from smallest to largest
  BS_mean = np.mean(BS_set)
  return ( BS_mean, BS_set[np.floor(bootNum*a1)-1], BS_set[np.ceil(bootNum*a2)-1] )

def LApermuPvalue(series1, series2, series3, precision, LA_score, fTransform, zNormalize):
  PP_set = np.zeros(precision, dtype='float')
  X = zNormalize(fTransform(series1))
  Y = zNormalize(fTransform(series2))
  Z = np.ma.array(series3)                                               #use = only assigns reference, must use a constructor
  for i in range(0, precision):
    np.random.shuffle(Z.T)
    PP_set[i] = compcore.calc_LA(X, Y, zNormalize(fTransform(Z)))
  if LA_score >= 0:
    P_two_tail = np.sum(np.abs(PP_set) >= LA_score)/float(precision)
  else:
    P_two_tail = np.sum(-np.abs(PP_set) <= LA_score)/float(precision)
  return P_two_tail
