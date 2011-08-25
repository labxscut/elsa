#lsalib.py -- Library of Local Similarity Analysis(LSA) Package
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

"""lsalib.py -- Library of Local Similarity Analysis(LSA) Package

  NOTE: numpy and scipy is required to use thie module
  NOTE: accepts input sequence table as delimited text file 
        with first row is the "Date" and other factor labels and
        first column is the time spot labels
"""

#import public resources
import csv, sys, random
import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.stats

#import lower level resource
try:
  #else run as installed
  from lsa import compcore
except ImportError:
  #try for debug
  import compcore

"""Synopsis:
    #Xf = simpleAverage( X )												#F: transform function
    #Xf = sdAverage( X )													#
    #Xz = scoreNormalize( X )												#N: normalize function
    #S_max = singleLSA( X, Y, D )											#D: delayLimit	
    #(S_low, S_up) = bootstrapCI( X, Y, D, alpha, B )						#B: bootNum, alpha: interval size
    #p = permuPvalue( X, Y, S_max, P )										#P: permuNum
    #q = storeyQvalue( p_set )												#p_set: set of p-values
    #lsaTable = applyAnalysis( rawData, delayLimit=D, ftransform=F, znormalize=N, bootCI=alpha, bootNum=B, permuNum=P )  
"""

###############################
# applyAnalsysi
# return a list-of-list, where the i-th list is:
#	[ f1, f2, L-score, L_low, L_up, X_start, Y_start, length, P/N, P-value, PCC,  PCC P-val,  Qvalue  ]
#	[ 1,  2,  3,       4,     5,    6,       7,       8,      9,   10,      11,   12,         13      ]
############################### 
        
def singleLSA(series1, series2, delayLimit, fTransform, zNormalize, keepTrace=True):
  """	do local simularity alignment 
		
		Args:
			series1(np.array): 	sequence data of Seq X
			series2(np.array): 	sequence data of Seq Y
			delayLimit(int): 	maximum time unit of delayed response allowed	
			fTransform(func):	replicate summarizing function

		Return:
			one single LSA result
      new implemented similarity alignment using external C++ routine in compcore
    
  """
  
  #print "f1=", fTransform(series1)
  #print "f2=", fTransform(series2)
  lsad=compcore.LSA_Data(delayLimit, zNormalize(fTransform(series1)), zNormalize(fTransform(series2)))
  lsar=compcore.DP_lsa(lsad, keepTrace)
  return lsar
	
def sample_wr(population, k):
  """ Chooses k random elements (with replacement) from a population
  """

  n = len(population)
  _random, _int = random.random, int  # speed hack 
  result = [None] * k
  for i in xrange(k):
    j = _int(_random() * n)
    result[i] = population[j]
  return np.array(result)

def bootstrapCI(series1, series2, Smax, delayLimit, bootCI, bootNum, fTransform, zNormalize, debug=0):
  """	do bootstrap CI estimation

		Args:
			series1(np.array): 	sequence data of Seq X
			series2(np.array): 	sequence data of Seq Y
			delayLimit(int): 	maximum time unit of delayed response allowed	
      bootCI(float):  confidence interval size
      bootNum(int): number of bootstraps
			fTransform(func):	replicate summarizing function

    Return:
      Confidence Interval

	"""

  ### no feasible, skipping bootstraping
  if series1.shape[0] == 1:
    return (Smax, Smax, Smax)

  ###print "------Bootstrapping------"
  lsad = compcore.LSA_Data()
  BS_set = np.zeros(bootNum, dtype='float')
  for i in range(0, bootNum):
    Xb = np.array([ sample_wr(series1[:,j], series1.shape[0]) for j in xrange(0,series1.shape[1]) ]).T
    Yb = np.array([ sample_wr(series2[:,j], series2.shape[0]) for j in xrange(0,series2.shape[1]) ]).T
    lsad.assign( delayLimit, zNormalize(fTransform(Xb)), zNormalize(fTransform(Yb)) )
    BS_set[i] = compcore.DP_lsa(lsad, False).score
  BS_set.sort()                                 #from smallest to largest
  BS_mean = np.mean(BS_set)
  #print np.histogram(BS_set, bins=10)
  a1 = (1-bootCI)/2.0
  a2 = bootCI+(1-bootCI)/2.0
  #bias correction steps
  if debug in [1, 3]:
    # correct Smax
    Smax = 2*Smax - BS_mean
  if debug in [2, 3]:
    # correct CI
    z0 = sp.stats.distributions.norm.ppf(np.sum(BS_set <= Smax)/float(bootNum))
    a1 = sp.stats.distributions.norm.cdf(2*z0+sp.stats.distributions.norm.ppf(a1))
    a2 = sp.stats.distributions.norm.cdf(2*z0+sp.stats.distributions.norm.ppf(a2))
    #print "z0=", z0, "a1=", a1, "a2=", a2
  return ( BS_mean, BS_set[np.floor(bootNum*a1)-1], BS_set[np.ceil(bootNum*a2)-1] )
	
def permuPvalue(series1, series2, delayLimit, permuNum, Smax, fTransform, zNormalize):
  """ do permutation Test

		Args:
			series1(np.array): 	sequence data of Seq X
			series2(np.array): 	sequence data of Seq Y
			delayLimit(int): 	maximum time unit of delayed response allowed	
      permuNum(int): number of permutations
      Smax(int): maximum LSA
			fTransform(func):	replicate summarizing function

    Return:
      p-value

	"""
	
  ###print "-------permutation------"
  lsad = compcore.LSA_Data()
  PP_set = np.zeros(permuNum, dtype='float')
  Xz = zNormalize(fTransform(series1))
  Y = np.array(series2)                                               #use = only assigns reference, must use a constructor
  for i in xrange(0, permuNum):
    np.random.shuffle(Y.T)
    lsad.assign( delayLimit, Xz, zNormalize(fTransform(Y)) )
    PP_set[i] = np.abs(compcore.DP_lsa(lsad, False).score)
  #PP_set[permuNum]=Smax                                               #the original test shall not be considerred
  #print "PP_set", PP_set, PP_set >= Smax, np.sum(PP_set>=Smax), np.float(permuNum)
  return np.sum(PP_set >= Smax)/np.float(permuNum)

def storeyQvalue(pvalues, lam=np.arange(0,0.9,0.05), method='smoother', robust=False, smooth_df=3):
  """ do Q-value calculation

    Args:
      pvalues(np.array):  a set of p-values
      lam(np.array):  tentative lambda data
      method(str):  calculating method, currently only support 'smoother'
      robust(bool): use robust static or not, default not
      smooth_df(int): order of spline function

    Returns:
      a set of qvalues
  """
  
  pi_set = np.zeros(len(lam), dtype='float')
  for i in xrange(0, len(lam)):
    pi_set[i] = np.mean(pvalues>=lam[i])/(1-lam[i])

  if method=='smoother':
    spline_fit = sp.interpolate.interp1d(lam, pi_set, kind=smooth_df)
    pi_0 = spline_fit(np.max(lam))
    pi_0 = np.max( [np.min( [np.min(pi_0), 1]), 0] ) #0<=pi_0<=1
    if pi_0 == 0:
      print >>sys.stderr, "Warning: smoother method not working, fall back to bootstrap"
      method='bootstrap'

  if method=='bootstrap':                            #bootstrap
    pi_min = np.min(pi_set)
    mse = np.zeros((100, len(lam)), dtype='float')
    pi_set_boot = np.zeros((100, len(lam)), dtype='float')
    for j in xrange(0, 100):
      p_boot = sample_wr(pvalues, len(pvalues))
      for i in xrange(0, len(lam)):
        pi_set_boot[j][i] = np.mean(p_boot>=lam[i])/(1-lam[i]) 
      mse[j] = (pi_set_boot[j]-pi_min)**2
    pi_0 = np.min(pi_set_boot[mse==np.min(mse)])
    pi_0 = np.max( [np.min( [np.min(pi_0), 1]), 0] ) #0<=pi_0<=1
    if pi_0 == 0:
      print >>sys.stderr, "Warning: bootstrap method not working, cannot estimate qvalues"
      return np.array( [np.nan] * len(pvalues) )

  #print "pvalues=", pvalues
  #print "pi_0=", pi_0
  p_argsort = np.argsort(pvalues)
  #print "argsort of p=", pvalues[p_argsort]
  p_ranks = tied_rank(pvalues)
  #print "tied rank of p=", p_ranks
  #print "lam,pi_set,pi_0:", lam, pi_set, pi_0
  #print "pi_0, p_ranks, pvalues, len(pvalues)", pi_0, p_ranks, pvalues, len(pvalues)
  if robust:
    qvalues = pi_0*len(pvalues)*pvalues/(p_ranks*(1-np.power((1-pvalues),len(pvalues))))
  else:
    qvalues = pi_0*len(pvalues)*pvalues/p_ranks  
  #print "qvalues=", qvalues
  qvalues[p_argsort[len(qvalues)-1]] = np.min( [qvalues[p_argsort[len(qvalues)-1]], 1] ) # argsort in asscending order
  for i in reversed(range(0,len(qvalues)-1)): #don't know why it is so complicated here, why not just use qvalues; to enssure desencing order!!!
    qvalues[p_argsort[i]] = np.min( [qvalues[p_argsort[i]], qvalues[p_argsort[i+1]], 1] )
  return qvalues

def simpleAverage(tseries):
  """ simple averaging 

    Args:
      tseries(np.array):  one time series with replicates, each row is a replicate

    Reterns:
      one row with replicates averaged

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
  """

  return np.average(np.nan_to_num(tseries), axis=0)

def sdAverage(tseries):
  """	SD weighted averaging 

    Args:
      tseries(np.array):  one time series with replicates, each row is a replicate

    Reterns:
      one row with replicates SD weighted averaged

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
  """
  sd = np.std(np.nan_to_num(tseries),axis=0,ddof=1)
  for v in sd:
    if v == 0:
      return np.average(np.nan_to_num(tseries), axis=0)                  #sd = 0, fall back to simpleAverage
  Xf = ((np.average(np.nan_to_num(tseries), axis=0)/sd)/np.sum(1/sd))/sd   #sd-weighted sample
  #return (Xf - np.average(Xf))/(np.sqrt(Xf.shape)*np.std(Xf))  #rescale and centralized
  return Xf

def simpleMedian(tseries):
  """ simple median

    Args:
      tseries(np.array):  one time series with replicates, each row is a replicate

    Reterns:
      one row with replicates summarized by median

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
  """

  return np.median(np.nan_to_num(tseries), axis=0)

def madMedian(tseries):
  """	MAD weighted averaging 

    Args:
      tseries(np.array):  one time series with replicates, each row is a replicate

    Reterns:
      one row with replicates summarized by MAD weighted median

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
  """
  Xf = np.nan_to_num(tseries)
  mad = np.median( np.abs(Xf - np.median(Xf, axis=0)), axis=0 )
  for v in mad:
    if v == 0:
      return np.median(np.nan_to_num(tseries), axis=0)                  #mad = 0, fall back to simpleMedian
  Xf = ((np.median(Xf, axis=0)/mad)/np.sum(1/mad))/mad                    #mad-weighted sample
  return Xf

def tied_rank(values):
  """ rank pvalues with ties, tie is ranked as the largest rank

    Args:
      values(np.array): a numeric array

    Returns:
      one vector of asscendant ranks from 1
       ties are kept and label by largest rank 
  """

  v_num = {}  #pvalues and counts
  v_cum = {} #cumulated rank, take largest for tie
  for v in values:
    v_num[v] = v_num.get(v,0) + 1
  suvs = v_num.keys()  #sorted unique pvaludes
  suvs.sort() 
  c = 0
  for v in suvs:
    c += v_num[v]
    v_cum[v] = c
  #print suvs
  #print v_cum
  #print [ v_cum[values[i]] for i in xrange(0, len(values)) ]
  return np.array( [ v_cum[values[i]] for i in xrange(0, len(values)) ], dtype='float' )


def	wholeNormalize(tseries):
  """	whole normalizing

    Args:
      tseries(np.array):  time series matrix with replicates

    Returns:
      wholely score normalized tseries
  """
  Xz = tseries              #make a copy
  #print "before normal, Xz=", Xz
  shape = Xz.shape          #save shape
  Xz = Xz.ravel()                #flatten
  #print Xz.ravel()
  nanb = np.isnan(Xz)       #find nans
  nans = (np.nonzero(nanb)[0]).tolist()   #save nans
  Xz = Xz[-nanb]            #cleaned, 
  ranks = tied_rank(Xz)     #rank na
  Xz = sp.stats.distributions.norm.ppf( ranks/(len(ranks)+1) )
  for nan in nans:
    Xz = np.insert(Xz, nan, 0.)   #insert zeros
  #print Xz
  #print Xz.shape, shape
  Xz.shape = shape
  #print "after normal, Xz=", Xz
  return Xz

def scoreNormalize(tseries):
  """ score normalizing

    Args:
      tseries(np.array): 1-d time series
    
    Returns:
      score normalized time series
  """
  ranks = tied_rank(tseries)
  return sp.stats.distributions.norm.ppf( ranks/(len(ranks)+1) )

def noneNormalize(tseries):
  """ no normalizaing

    Args:
      tseries(np.array):  time series matrix

    Returns:
      non normalized tseries
  """

  return np.nan_to_num(tseries) 


def fillMissing(tseries, method):
  """ fill missing data

    Args:
      tseries(np.array):  one time series with no replicates
      method(str):  filling method, choices ['none', 'zero', 'linear', 'slinear', 'nearest', 'quadratic', 'cubic'] 

    Reterns:
      tseries with missing data filled
  """
  
  if method == 'none':
    return tseries #return with nan's
  else:
    y = tseries[np.logical_not(np.isnan(tseries))]
    x = np.array(range(0, len(tseries)))[np.logical_not(np.isnan(tseries))]
    try:
      spline_fit = sp.interpolate.interp1d( x, y, method )
    except:
      print >>sys.stderr, "cannot fill missing values using ", method, "method, fall back to none" 
      return np.nan_to_num(tseries)
    yy = np.zeros( len(tseries), dtype='float' )
    for i in range(0, len(tseries)):
      try:
        yy[i] = spline_fit(i)
      except ValueError:
        yy[i] = 0
    return yy
    
def applyAnalysis(cleanData, delayLimit=3, bootCI=.95, bootNum=1000, permuNum=1000, fTransform=simpleAverage, zNormalize=scoreNormalize):
  """ calculate pairwise LS scores and p-values

    	Args:
    		cleanData(np.array): 	numpy data array with correct format factor_num x timespot_num x replicte_num, possibly nans
    		delayLimit(int): 	maximum time unit of delayed response allowed
     		bootCI(float): 		bootstrap confidence interval size, 0 to 1
    		bootNum(int): 		bootstrap number
    		permuNum(int): 		number of permutations
    		ftransform(func): 	summarizing function for replicated data
    		znormalize(func): 	normalizing function for ftransformed data
    		
    	Returns:
    		A LSA table.
    		each row in the table is a list in following format:
    		[ Seq X's Idx, Seq Y's Idx, LS Score, CI_low, CI_high, X's Start Position, 
        	Y's Start Position, Alignment Length, X delay to Y,
        	P-value, Pearson' Correlation, P-value of PCC, Q-value ]
        	
  """	

  factorNum = cleanData.shape[0]
  repNum = cleanData.shape[1]
  spotNum = cleanData.shape[2]
  lsaTable = [None]*(factorNum*(factorNum-1)/2)
  pvalues = np.zeros(factorNum*(factorNum-1)/2, dtype='float')
  pccpvalues = np.zeros(factorNum*(factorNum-1)/2, dtype='float')

  #print factorNum, repNum, spotNum, lsaTable, pvalues

  ti = 0
  for i in xrange(0, factorNum-1):
    Xz = cleanData[i]
    for j in xrange(i+1, factorNum):
      #print "normalizing..."
      Yz = cleanData[j]
      #print "lsa computing..."
      #print "Xz=", Xz
      #print "Yz=", Yz
      LSA_result = singleLSA(Xz, Yz, delayLimit, fTransform, zNormalize, True)                          # do LSA computation
      Smax = LSA_result.score                                                               # get Smax
      #print "bootstrap computing..."
      (Smax, Sl, Su) = bootstrapCI(Xz, Yz, Smax, delayLimit, bootCI, bootNum, fTransform, zNormalize)           # do Bootstrap CI
      #print "permutation test..."
      permuP = permuPvalue(Xz, Yz, delayLimit, permuNum, np.abs(Smax), fTransform, zNormalize)          # do Permutation Test
      pvalues[ti] = permuP
      Al = len(LSA_result.trace)
      (Xs, Ys, Al) = (LSA_result.trace[Al-1][0], LSA_result.trace[Al-1][1], len(LSA_result.trace))
      #print "PPC..." 
      #print np.mean(Xz, axis=0), np.mean(Yz, axis=0)
      (PCC, P_PCC) = sp.stats.pearsonr(np.mean(np.nan_to_num(Xz), axis=0), np.mean(np.nan_to_num(Yz), axis=0))
      # need +epsilon to avoid all zeros
      pccpvalues[ti] = P_PCC
      lsaTable[ti] = [i, j, Smax, Sl, Su, Xs, Ys, Al, Xs-Ys, permuP, PCC, P_PCC]
      ti += 1
      #print "finalizing..."

  #print "qvalues ..."
  #print "pvalues", pvalues
  qvalues = storeyQvalue( pvalues )
  pccqvalues = storeyQvalue( pccpvalues )
  for i in xrange(0, len(qvalues)):
    lsaTable[i].append( qvalues[i] )
    lsaTable[i].append( pccqvalues[i] )

  return lsaTable

def test():
  """ self test script
  """
  print >>sys.stderr, "###testing###"
  test_data = np.array( [ [ [ 3, 2, 3, 4], [1, 6, 2, 3], [6, 4, 6, 8] ], [ [np.nan, 2, np.nan, 3], [4, 5, 3, np.nan], [1, 1, 1, 1] ] ], dtype='float' )
  test_fN = test_data.shape[0]
  test_rN = test_data.shape[1]
  test_tN = test_data.shape[2]
  print >>sys.stderr, "fN, tN, rN", test_fN, test_tN, test_rN
  print >>sys.stderr, "test_data", test_data
  print >>sys.stderr, "---fillMissing---"
  print >>sys.stderr, "none:", fillMissing(test_data[1][0], 'none')
  print >>sys.stderr, "zero:", fillMissing(test_data[1][0], 'zero')
  print >>sys.stderr, "linear:", fillMissing(test_data[1][0], 'linear')
  print >>sys.stderr, "quadratic:", fillMissing(test_data[1][0], 'quadratic')
  print >>sys.stderr, "cubic:", fillMissing(test_data[1][0], 'cubic')
  print >>sys.stderr, "slinear:", fillMissing(test_data[1][0], 'slinear')
  print >>sys.stderr, "nearest:", fillMissing(test_data[1][0], 'nearest')
  print >>sys.stderr, "###use zero as clean data###"
  clean_data = test_data
  for i in xrange(0, test_fN):
    for j in xrange(0, test_rN):
      clean_data[i][j] = fillMissing(test_data[i][j], 'none')
  print >>sys.stderr, "clean_data", clean_data
  print >>sys.stderr, "---tied-rank---"
  print >>sys.stderr, tied_rank(clean_data[0][0])
  print >>sys.stderr, "---scoreNormalize---"
  print >>sys.stderr, scoreNormalize(clean_data[0][0])
  print >>sys.stderr, scoreNormalize(clean_data[0][1])
  print >>sys.stderr, "---wholeNormalize---"
  print >>sys.stderr, wholeNormalize(clean_data[0])
  print >>sys.stderr, wholeNormalize(clean_data[0])
  print >>sys.stderr, "---simpleAverage---" 
  print >>sys.stderr, simpleAverage(clean_data[1])
  print >>sys.stderr, "---sdAverage---"
  print >>sys.stderr, sdAverage(clean_data[1])
  print >>sys.stderr, "---simpleMedian---" 
  print >>sys.stderr, simpleMedian(clean_data[1])
  print >>sys.stderr, "---madMedian---"
  print >>sys.stderr, madMedian(clean_data[1])
  print >>sys.stderr, "---storeyQvalue---"
  pvalues = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.02, 0.03, 0.04, 0.03, 0.03], dtype='float')
  print >>sys.stderr, "pvalues:", pvalues 
  print >>sys.stderr, storeyQvalue(pvalues)
  print >>sys.stderr, "---singleLSA---"
  lsar = singleLSA(clean_data[0], clean_data[1], delayLimit=1, fTransform=simpleAverage, zNormalize=scoreNormalize, keepTrace=True)
  print >>sys.stderr, lsar.score
  print >>sys.stderr, "---bootstrapCI---"
  print >>sys.stderr, bootstrapCI(clean_data[0], clean_data[1], lsar.score, 1, .95, 10, simpleAverage, scoreNormalize)
  print >>sys.stderr, "---permuPvalue---"
  print >>sys.stderr, "p-value:", permuPvalue(clean_data[1], clean_data[0], 1, 10, np.abs(lsar.score), simpleAverage, scoreNormalize)
  print >>sys.stderr, "---PCC---"
  (nPCC, nP_PCC) = sp.stats.pearsonr(np.mean(np.nan_to_num(test_data[0]), axis=0), np.mean(np.nan_to_num(test_data[1]), axis=0))
  oPCC = sp.corrcoef( np.mean(np.nan_to_num(test_data[0]),axis=0), np.mean(np.nan_to_num(test_data[1]),axis=0) )[0,1]
  otcdf = sp.stats.distributions.t.cdf(oPCC*np.sqrt((test_tN-2)/np.float(1.000000001-oPCC**2)), (test_tN-2))
  oP_PCC = .5 + np.sign(oPCC)*(.5 - otcdf) #addhoc for perfect correlation
  print >>sys.stderr, "nPCC", "nP_PCC", "oPCC", "oP_PCC", "otcdf"
  print >>sys.stderr, nPCC, nP_PCC, oPCC, oP_PCC, otcdf
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, 1, .95, 10, 10, simpleAverage, scoreNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, 1, .95, 10, 10, sdAverage, scoreNormalize)
  

if __name__=="__main__":
  print "hello world!"
  test()
  exit(0)


#############  FUNCTIONS NO IN USE AFTER THIS LINE ####################

