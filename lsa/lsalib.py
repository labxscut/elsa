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
#import numpy.ma as np.ma
import scipy as sp
import scipy.interpolate
import scipy.stats

#import lower level resource
try:
  #else run as installed
  from lsa import compcore
  np.seterr(all='ignore')                             #ignore RuntimeWarning of abs in installed mode
except ImportError:
  #try for debug
  import compcore
  np.seterr(all='warn')


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
# applyAnalsys
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
  result = np.array([np.nan] * k)
  for i in xrange(k):
    j = _int(_random() * n)
    if type(population) == np.ma.MaskedArray:
      if population.mask[j]:
        result[i] = np.nan
      else:
        result[i] = population[j]
    else:
      result[i] = population[j]
  if type(population) == np.ma.MaskedArray:
    result = np.ma.masked_invalid(result)
  return result

def ma_median(ts, axis=0):
  ns = np.ma.median(ts, axis=axis)
  if type(ns.mask) == np.bool_:       #fix broken ma.median, mask=False instead of [False, ...] for all mask
    ns.mask = [ns.mask] * ns.shape[axis]
  return ns

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
    Xb = np.ma.array([ sample_wr(series1[:,j], series1.shape[0]) for j in xrange(0,series1.shape[1]) ]).T
    Yb = np.ma.array([ sample_wr(series2[:,j], series2.shape[0]) for j in xrange(0,series2.shape[1]) ]).T
    #print "Xb=", Xb
    #print "Yb=", Yb
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

def theoPvalue(xx, timespots, D, Kcut):
  A = (1/xx)
  B = 2*D+1
  P = 1
  for k in xrange(1,Kcut+1):
    B = (2*k-1)**2
    P = P-(8**B)*((A+pipi_inv*(1/C))*exp(-C*pipi/(2*xx))**B)
  #if P <= 0:
  #  P = 0.
  return P
	
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
  Y = np.ma.array(series2)                                               #use = only assigns reference, must use a constructor
  for i in xrange(0, permuNum):
    np.random.shuffle(Y.T)
    lsad.assign( delayLimit, Xz, zNormalize(fTransform(Y)) )
    PP_set[i] = compcore.DP_lsa(lsad, False).score
  #PP_set[permuNum]=Smax                                               #the original test shall not be considerred
  #print "PP_set", PP_set, PP_set >= Smax, np.sum(PP_set>=Smax), np.float(permuNum)
  if Smax >= 0:
    P_one_tail = np.sum(PP_set >= Smax)/np.float(permuNum)
  else:
    P_one_tail = np.sum(PP_set <= Smax)/np.float(permuNum)
  return P_one_tail

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
  #print "data=", tseries
  #print "simpleAverage=", np.ma.average(tseries, axis=0)
  return np.ma.average(tseries, axis=0)

def sdAverage(tseries):
  """	SD weighted averaging 

    Args:
      tseries(np.array):  one time series with replicates, each row is a replicate

    Reterns:
      one row with replicates SD weighted averaged

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
  """
  sd = np.ma.std(tseries,axis=0,ddof=1)
  for v in sd:
    if v == 0.:
      return np.ma.average(tseries, axis=0)                       #sd = 0, fall back to simpleAverage
  #print "tseries=", tseries
  #print "sd=", sd
  Xf = ((np.ma.average(tseries, axis=0)/sd)/np.ma.sum(1/sd))/sd   #sd-weighted sample
  #return (Xf - np.average(Xf))/(np.sqrt(Xf.shape)*np.std(Xf))  #rescale and centralized
  #print "sdAverage=", Xf
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
  
  Xf = ma_median(tseries, axis=0)
  #print "simpleMedian=", Xf, Xf.mask, type(Xf.mask) #ma.median broken
  return Xf

def madMedian(tseries):
  """	MAD weighted averaging 

    Args:
      tseries(np.array):  one time series with replicates, each row is a replicate

    Reterns:
      one row with replicates summarized by MAD weighted median

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
  """
  Xf = tseries
  mad = ma_median( np.ma.abs(Xf - ma_median(Xf, axis=0)), axis=0 ) #current throw a warning 
  for v in mad:
    if v == 0:
      return ma_median(tseries, axis=0)                  #mad = 0, fall back to simpleMedian
  Xf = ((ma_median(Xf, axis=0)/mad)/np.ma.sum(1/mad))/mad                    #mad-weighted sample

  #print "medMedian=", Xf, Xf.mask, type(Xf.mask)
  return Xf

def tied_rank(values):
  """ rank pvalues with ties, tie is ranked as the largest rank, now allow na's non-ranked

    Args:
      values(np.array): a numeric array

    Returns:
      one vector of asscendant ranks from 1
       ties are kept and label by largest rank 
  """

  V = values
  #print "nonzero(V.mask)=", np.nonzero(V.mask)
  #print "V=", V
  if type(V) == np.ma.MaskedArray:
    #print "V.mask=", V.mask
    nans = (np.nonzero(V.mask)[0]).tolist()         #record indecies
    V = V[-V.mask]                               #remove nan's
    #print "V=", V

  v_num = {}  #values and counts
  v_cum = {}  #cumulated rank, take largest for tie
  for v in V:
    v_num[v] = v_num.get(v,0) + 1
  suvs = v_num.keys()  #sorted unique pvaludes
  suvs.sort() 
  c = 0
  for v in suvs:
    c += v_num[v]
    v_cum[v] = c
  #print "sorted unique values=", suvs
  #print "cumulative number=", v_cum
  #print [ v_cum[values[i]] for i in xrange(0, len(values)) ]
  sV = np.array( [ v_cum[V[i]] for i in xrange(0, len(V)) ], dtype='float' ) #sorted V
  #print "sV=", sV
  #print "nans=", nans
  
  if type(V) == np.ma.MaskedArray:
    #mV = sV
    for nan in nans:
      sV = np.insert(sV, nan, np.nan)   #insert nan to original place, need return masked? yes
    sV = np.ma.masked_invalid(sV)
    #print "sV=", sV

  return sV 

#def	wholeNormalize(tseries):
#  """	whole normalizing
#
#    Args:
#      tseries(np.array):  time series matrix with replicates
#
#    Returns:
#      wholely score normalized tseries
#  """
#  Xz = tseries              #make a copy
#  #print "before normal, Xz=", Xz
#  shape = Xz.shape          #save shape
#  Xz = Xz.ravel()                #flatten
#  #print Xz.ravel()
#  ranks = tied_rank(Xz)     #rank na
#
#  Xz = sp.stats.distributions.norm.ppf( ranks/(len(ranks)+1) )
#  #print Xz
#  #print Xz.shape, shape
#  Xz.shape = shape
#  #print "after normal, Xz=", Xz
#  return Xz

def scoreNormalize(tseries):
  """ score normalizing

    Args:
      tseries(np.array): 1-d time series
    
    Returns:
      score normalized time series
  """
  ranks = tied_rank(tseries)
  nt = sp.stats.distributions.norm.ppf( ranks/(len(ranks)+1) )
  #print "nt=", nt
  nt = np.nan_to_num(nt)  #filling zeros to nan, shall be no na's from here 
  return nt

def noneNormalize(tseries):
  """ no normalizaing

    Args:
      tseries(np.array):  time series matrix

    Returns:
      non normalized tseries
  """

  return np.nan_to_num(tseries)  #filling zeros to nan

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
    
def applyAnalysis(firstData, secondData, onDiag=True, delayLimit=3, bootCI=.95, bootNum=1000, permuNum=1000, fTransform=simpleAverage, zNormalize=scoreNormalize):
  """ calculate pairwise LS scores and p-values

    	Args:
    		firstData(np.array): 	numpy data array with correct format factor_num x timespot_num x replicte_num, possibly nans
    		secondData(np.array): 	numpy data array with correct format factor_num x timespot_num x replicte_num, possibly nans
                noDiag(bool):           no results for diagnol comparisions
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

  pipi = math.pi**2 # pi^2
  pipi_inv = 1/pipi

  firstFactorNum = firstData.shape[0]
  firstRepNum = firstData.shape[1]
  firstSpotNum = firstData.shape[2]
  secondFactorNum = secondData.shape[0]
  secondRepNum = secondData.shape[1]
  secondSpotNum = secondData.shape[2]
  #for now let's assume it is square
  #assert secondFactorNum == firstFactorNum and secondRepNum == firstRepNum and secondSpotNum == firstSpotNum 
  if onDiag:  # if assigned jobs are on the diagnal
    assert firstFactorNum == secondFactorNum
    pairwiseNum = firstFactorNum*(firstFactorNum-1)/2
  else:
    pairwiseNum = firstFactorNum*scecondFactorNum
  lsaTable = [None]*pairwiseNum
  pvalues = np.zeros(pairwiseNum, dtype='float')
  pccpvalues = np.zeros(pairwiseNum, dtype='float')
  #print factorNum, repNum, spotNum, lsaTable, pvalues

  ti = 0
  for i in xrange(0, firstFactorNum):
    Xz = np.ma.masked_invalid(firstData[i]) #need to convert to masked array with na's
    for j in xrange(0, secondFactorNum):
      if onDiag and i<=j:
        continue   #only care lower triangle entries
      #print "normalizing..."
      Yz = np.ma.masked_invalid(secondData[j]) # need to convert to masked array with na's
      timespots = Yz.shape[1]
      #print "lsa computing..."
      #print "Xz=", Xz
      #print "Yz=", Yz
      LSA_result = singleLSA(Xz, Yz, delayLimit, fTransform, zNormalize, True)                          # do LSA computation
      Smax = LSA_result.score                                                               # get Smax
      Al = len(LSA_result.trace)
      (Xs, Ys, Al) = (LSA_result.trace[Al-1][0], LSA_result.trace[Al-1][1], len(LSA_result.trace))
      #print "bootstrap computing..."
      if bootNum > 0: #do BS
        (Smax, Sl, Su) = bootstrapCI(Xz, Yz, Smax, delayLimit, bootCI, bootNum, fTransform, zNormalize)           # do Bootstrap CI
      else: #skip BS
        (Smax, Sl, Su) = (Smax, Smax, Smax)
      #print "permutation test..."
      if permuNum >= 0:
        lsaP = permuPvalue(Xz, Yz, delayLimit, permuNum, np.abs(Smax), fTransform, zNormalize)          # do Permutation Test
      else:
        xx = (np.abs(Smax)*np.sqrt(timespots))**2
        pipi_over_xx = pipi/xx
        D = np.abs(Xs-Ys)
        alpha = 1/permuNum/10
        Kcut = np.ceil( .5 - np.log( alpha/(8**(2*D+1))*xx*(1-exp(-pipi_over_xx))/2 )/pipi_over_xx )
        print Kcut
        lsaP = theoPvalue(xx, D, Al, Kcut)  #Rn/sqrt(n)=Smax*sqrt(n)
      pvalues[ti] = lsaP
      #print "PPC..." 
      #print np.mean(Xz, axis=0), np.mean(Yz, axis=0)
      (PCC, P_PCC) = sp.stats.pearsonr(np.ma.average(Xz, axis=0), np.ma.average(Yz, axis=0)) # two tailed p-value
      P_PCC = P_PCC/2   # one tailed p-value
      # need +epsilon to avoid all zeros
      pccpvalues[ti] = P_PCC
      lsaTable[ti] = [i, j, Smax, Sl, Su, Xs, Ys, Al, Xs-Ys, lsaP, PCC, P_PCC]
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
  np.seterr(all='warn')
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
  print >>sys.stderr, "---masked_array---"
  masked_data = np.ma.masked_invalid(clean_data)
  print >>sys.stderr, "masked_data", masked_data
  print >>sys.stderr, "---tied-rank---"
  print >>sys.stderr, "data:", masked_data[1][0]
  print >>sys.stderr, tied_rank(masked_data[1][0])
  print >>sys.stderr, "data:", masked_data[1][1]
  print >>sys.stderr, tied_rank(masked_data[1][1])
  #print >>sys.stderr, "---wholeNormalize---"
  #print >>sys.stderr, wholeNormalize(clean_data[0])
  #print >>sys.stderr, wholeNormalize(clean_data[0])
  print >>sys.stderr, "---simpleAverage---" 
  print >>sys.stderr, simpleAverage(masked_data[1])
  print >>sys.stderr, "---sdAverage---"
  print >>sys.stderr, sdAverage(masked_data[1])
  print >>sys.stderr, "---simpleMedian---" 
  print >>sys.stderr, simpleMedian(masked_data[1])
  print >>sys.stderr, "---madMedian---"
  print >>sys.stderr, madMedian(masked_data[1])
  print >>sys.stderr, "---scoreNormalize---"
  print >>sys.stderr, scoreNormalize(masked_data[1][0])
  print >>sys.stderr, scoreNormalize(masked_data[1][1])
  print >>sys.stderr, "---storeyQvalue---"
  pvalues = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.02, 0.03, 0.04, 0.03, 0.03], dtype='float')
  print >>sys.stderr, "pvalues:", pvalues 
  print >>sys.stderr, storeyQvalue(pvalues)
  print >>sys.stderr, "---singleLSA---"
  lsar = singleLSA(masked_data[0], masked_data[1], delayLimit=1, fTransform=simpleAverage, zNormalize=scoreNormalize, keepTrace=True)
  print >>sys.stderr, "lsar.score=", lsar.score
  print >>sys.stderr, "---bootstrapCI---"
  print >>sys.stderr, "Bset=", bootstrapCI(masked_data[0], masked_data[1], lsar.score, 1, .95, 2, simpleAverage, scoreNormalize)
  print >>sys.stderr, "---permuPvalue---"
  print >>sys.stderr, "P=", permuPvalue(masked_data[1], masked_data[0], 1, 2, np.abs(lsar.score), simpleAverage, scoreNormalize)
  print >>sys.stderr, "---PCC---"
  (nPCC, nP_PCC) = sp.stats.pearsonr(np.mean(np.nan_to_num(test_data[0]), axis=0), np.mean(np.nan_to_num(test_data[1]), axis=0))
  oPCC = sp.corrcoef( np.mean(np.nan_to_num(test_data[0]),axis=0), np.mean(np.nan_to_num(test_data[1]),axis=0) )[0,1]
  otcdf = sp.stats.distributions.t.cdf(oPCC*np.sqrt((test_tN-2)/np.float(1.000000001-oPCC**2)), (test_tN-2))
  oP_PCC = .5 + np.sign(oPCC)*(.5 - otcdf) #addhoc for perfect correlation
  print >>sys.stderr, "nPCC", "nP_PCC", "oPCC", "oP_PCC", "otcdf"
  print >>sys.stderr, nPCC, nP_PCC, oPCC, oP_PCC, otcdf
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, 1, .95, 10, 10, sdAverage, scoreNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, 1, .95, 10, 10, simpleAverage, scoreNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, 1, .95, 10, 10, simpleMedian, scoreNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, 1, .95, 10, 10, madMedian, scoreNormalize)

if __name__=="__main__":
  print "hello world!"
  test()
  exit(0)


#############  FUNCTIONS NO IN USE AFTER THIS LINE ####################

