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
    #p = permuPvalue( X, Y, S_max, P )										#P: pvalueMethod
    #q = storeyQvalue( p_set )												#p_set: set of p-values
    #lsaTable = applyAnalysis( rawData, delayLimit=D, ftransform=F, znormalize=N, bootCI=alpha, bootNum=B, pvalueMethod=P )  
"""

#global variable, stores calculated p-values.
P_table = dict()
my_decimal = 2    # preset x step size for P_table
pipi = np.pi**2 # pi^2
pipi_inv = 1/pipi

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

def readPvalue(S, timespots, x_decimal=2):
  # x = S*sqrt(timespots) => S=x/sqrt(timespots);
  # has to ceil the x value to avoid round to 0, which is not amenable to calculation
  #print "x=", int(np.around(S*np.sqrt(timespots)*(10**x_decimal)))
  xi = int(np.around(S*np.sqrt(timespots)*(10**x_decimal)))
  if xi in P_table:
    return P_table[xi]
  else:
    return 0.

def theoPvalue(timespots, Dmax, precision, x_decimal=2):   #let's produce 2 tail-ed p-value
  #print np.linspace(0,timespots,timespots/10**(-x_decimal)+1)
  for xi in xrange(0,timespots*(10**(x_decimal))+1): 
    if xi == 0:
      P_table[xi] = 1
      continue
    #print x
    x = xi/float(10**(x_decimal))
    xx = x**2
    pipi_over_xx = pipi/xx
    alpha = precision
    B = 2*Dmax+1
    #print pipi_over_xx
    #print alpha
    #print np.log(alpha*xx*(1-np.exp(-pipi_over_xx))/(8**B)/2)
    #print np.log(alpha*xx*(1-np.exp(-pipi_over_xx))/(8**B)/2) / pipi_over_xx
    Kcut = int(np.ceil( .5 - np.log( (alpha/(2**B-1))**(1/B) *xx*(1-np.exp(-pipi_over_xx))/8/2 )/pipi_over_xx ))
    #Kcut = 200
    A = 1/xx
    R = 0     # root of cdf

    for k in xrange(1,Kcut+1):
      C = (2*k-1)**2
      R = R + (A+pipi_inv/C)*np.exp(-C*pipi_over_xx/2)
      P_two_tail = 1 - (8**B)*(R**B)
      #if x==2.52:
      #  print "k=",k, "A=",A, "B=",B, "C=",C, "F=",(8**B)*(R**B), "dR=",(A+pipi_inv/C)*np.exp(-C*pipi_over_xx/2), "P=",P,"pipi_inv=",pipi_inv,"pipi_over_xx=",pipi_over_xx 
    #print xi, x, Kcut, P, P/2;

    if P_two_tail < precision:
      P_table[xi] = P_two_tail
      break
    else:
      P_table[xi] = P_two_tail #return one tailed probability

  return
	
def permuPvalue(series1, series2, delayLimit, pvalueMethod, Smax, fTransform, zNormalize):
  """ do permutation Test

    Args:
			series1(np.array): 	sequence data of Seq X
			series2(np.array): 	sequence data of Seq Y
			delayLimit(int): 	maximum time unit of delayed response allowed	
      pvalueMethod(int): number of permutations
      Smax(int): maximum LSA
			fTransform(func):	replicate summarizing function

    Return:
      p-value

	"""
	
  ###print "-------permutation------"
  lsad = compcore.LSA_Data()
  PP_set = np.zeros(pvalueMethod, dtype='float')
  Xz = zNormalize(fTransform(series1))
  Y = np.ma.array(series2)                                               #use = only assigns reference, must use a constructor
  for i in xrange(0, pvalueMethod):
    np.random.shuffle(Y.T)
    lsad.assign( delayLimit, Xz, zNormalize(fTransform(Y)) )
    PP_set[i] = compcore.DP_lsa(lsad, False).score
  #PP_set[pvalueMethod]=Smax                                               #the original test shall not be considerred
  #print "PP_set", PP_set, PP_set >= Smax, np.sum(PP_set>=Smax), np.float(pvalueMethod)
  if Smax >= 0:
    P_two_tail = np.sum(np.abs(PP_set) >= Smax)/np.float(pvalueMethod)
  else:
    P_two_tail = np.sum(np.abs(PP_set) <= Smax)/np.float(pvalueMethod)
  return P_two_tail

def storeyQvalue(pvalues, lam=np.arange(0,0.9,0.05), method='smoother', robust=False, smooth_df=3):
  """ do Q-value calculation

    Args:
      pvalues(np.array):  a set of p-values
      lam(np.array):  tentative lambda data
      method(str):  calculating method, currently only support 'smoother'
      robust(bool): use robust static or not, default not
      smooth_df(int): order of spline function

    Returns:
      qvalues(np.array): a set of qvalues
  """
  
  #assume pvalues is an array with possible nans
  #get a list of non nan pvalues, do the same procedure, putback to original list
  mpvalues = np.ma.masked_invalid(pvalues)
  rpvalues = mpvalues[~mpvalues.mask]
  p_num = len(pvalues)
  rp_num = len(rpvalues)
  #print "p_num=", p_num, "rp_num=", rp_num
  #print "rpvalues=", rpvalues

  if rp_num <= 1:
    #print >>sys.stderr, "WARN: not enough number of pvalues for q-value evaluation! nan will be filled!"
    return np.array( [np.nan] * p_num, dtype='float')

  pi_set = np.zeros(len(lam), dtype='float')
  for i in xrange(0, len(lam)):
    pi_set[i] = np.mean(rpvalues>=lam[i])/(1-lam[i])

  #print "pi_set=", pi_set
  if method=='smoother':
    spline_fit = sp.interpolate.interp1d(lam, pi_set, kind=smooth_df)
    pi_0 = spline_fit(np.max(lam))
    #print "pi_0=", pi_0
    pi_0 = np.max( [np.min( [np.min(pi_0), 1]), 0] ) #0<=pi_0<=1
    #print "pi_0=", pi_0
    if pi_0 == 0:
      #print >>sys.stderr, "WARN: smoother method not working, fall back to bootstrap"
      method='bootstrap'

  if method=='bootstrap':                            #bootstrap
    pi_min = np.min(pi_set)
    mse = np.zeros((100, len(lam)), dtype='float')
    pi_set_boot = np.zeros((100, len(lam)), dtype='float')
    for j in xrange(0, 100):
      p_boot = sample_wr(rpvalues, rp_num)
      for i in xrange(0, len(lam)):
        pi_set_boot[j][i] = np.mean(p_boot>=lam[i])/(1-lam[i]) 
      mse[j] = (pi_set_boot[j]-pi_min)**2
    min_mse_j = np.argmin(mse)
    pi_0 = np.min(pi_set_boot[min_mse_j])
    #print "pi_0=", pi_0
    pi_0 = np.max([np.min( [np.min(pi_0), 1]), 0]) #0<=pi_0<=1
    #print "pi_0=", pi_0
    if pi_0 == 0:
      #print >>sys.stderr, "WARN: bootstrap method not working, cannot estimate qvalues"
      return np.array( [np.nan] * p_num, dtype='float' )

  #print "pi_0=", pi_0
  rp_argsort = np.argsort(rpvalues)                     #np.nan will be sorted as maximum (the largest rank)
  #print "argsort of rps=", rpvalues[rp_argsort]
  rp_ranks = tied_rank(rpvalues)
  #print "tied rank of rps=", rp_ranks
  #print "lam,pi_set,pi_0:", lam, pi_set, pi_0
  #print "pi_0, p_ranks, pvalues, len(pvalues)", pi_0, p_ranks, pvalues, len(pvalues)
  if robust:
    rqvalues = pi_0*rp_num*rpvalues*(1/(rp_ranks*(1-np.power((1-rpvalues),rp_num))))
  else:
    rqvalues = pi_0*rp_num*rpvalues*(1/rp_ranks) 
  #print "rqs=", rqvalues
  rqvalues[rp_argsort[rp_num-1]] = np.min( [rqvalues[rp_argsort[rp_num-1]], 1] ) # argsort in asscending order
  for i in reversed(range(0,rp_num-1)): #don't know why it is so complicated here, why not just use qvalues; to enssure desencing order!!!
    rqvalues[rp_argsort[i]] = np.min( [rqvalues[rp_argsort[i]], rqvalues[rp_argsort[i+1]], 1] )

  qvalues=np.array([np.nan]*p_num)
  j=0
  for i in range(0, p_num):
    if not mpvalues.mask[i]:
      qvalues[i]=rqvalues[j]
      j += 1

  #print "qs=", qvalues
  return qvalues

def simpleAverage(tseries):
  """ simple averaging 

    Args:
      tseries(np.ma.array):  one 2d time series (masked array) with replicates, each row is a replicate

    Reterns:
      (1d np.ma.array) one row with replicates averaged

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
  """
  #print "data=", tseries
  #print "simpleAverage=", np.ma.average(tseries, axis=0)
  return np.ma.average(tseries, axis=0)

def sdAverage(tseries):
  """	SD weighted averaging 

    Args:
      tseries(np.ma.array):  one 2d time series (masked array) with replicates, each row is a replicate

    Reterns:
      (1d np.ma.array): one row with replicates SD weighted averaged

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
  """
  #print tseries
  try:
    sd = np.ma.std(tseries,axis=0,ddof=1)
  except FloatingPointError:
    return np.ma.average(tseries, axis=0)                       #sd = 0, fall back to simpleAverage
  if np.any(sd.mask) or (np.ma.sum(sd==0))>0:
    #print sd, sd.mask, sd==0
    #print np.ma.average(tseries, axis=0)
    return np.ma.average(tseries, axis=0)                       #sd = 0, fall back to simpleAverage
  #print "tseries=", tseries
  #print "sd=", sd
  #print "average=", np.ma.average(tseries, axis=0)
  #av_inv_sd = np.ma.average(tseries, axis=0)/sd
  #sum_inv_sd = np.ma.sum(1/sd)
  #print "sum(1/sd)=", np.ma.sum(1/sd)
  #av_inv_sum = av_inv_sd/sum_inv_sd
  #print "av/sum(1/sd)", av_inv_sum
  #print "inv_sd=", 1/sd
  #Xf = av_inv_sum*(1/sd)
  #print "Xf=", Xf
  Xf = np.ma.average(tseries, axis=0)*(1/sd)*(1/np.ma.sum(1/sd))*(1/sd)   #sd-weighted sample
  #print "Xf=", Xf
  #Xf = np.divide(np.divide(np.divide(np.ma.average(tseries, axis=0),sd),np.ma.sum(1/sd)),sd)   #sd-weighted sample
  #Xf = ((np.ma.average(tseries, axis=0)/sd)/np.ma.sum(1/sd))/sd   
  #wierd warning: ma/core.py:772: RuntimeWarning: underflow encountered in multiply return umath.absolute(a) * self.tolerance >= umath.absolute(b)
  #return (Xf - np.average(Xf))/(np.sqrt(Xf.shape)*np.std(Xf))  #rescale and centralized
  #print "sdAverage=", Xf
  return Xf

def simpleMedian(tseries):
  """ simple median

    Args:
      tseries(2d np.mp.array):  one time series with replicates, each row is a replicate

    Reterns:
      1d np.mp.array: one row with replicates summarized by median

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
  """
  
  Xf = ma_median(tseries, axis=0)
  #print "simpleMedian=", Xf, Xf.mask, type(Xf.mask) #ma.median broken
  return Xf

def madMedian(tseries):
  """	MAD weighted averaging 

    Args:
      tseries(2d np.ma.array):  one time series with replicates, each row is a replicate

    Reterns:
      1d np.ma.array: one row with replicates summarized by MAD weighted median

    Note:
      if nan in tseries, it is treated as zeros, this will happen if fTransform before zNormalize
  """
  Xf = tseries
  mad = ma_median( np.ma.abs(Xf - ma_median(Xf, axis=0)), axis=0 ) #current throw a warning 
  #print mad, mad.mask, mad==0
  if np.any(mad.mask) or (np.ma.sum(mad==0))>0:
    return ma_median(tseries, axis=0)                  #mad = 0, fall back to simpleMedian
  Xf = ma_median(Xf, axis=0)*(1/mad)*(1/np.ma.sum(1/mad))*(1/mad)                   #mad-weighted sample

  #print "medMedian=", Xf, Xf.mask, type(Xf.mask)
  return Xf

def tied_rank(values):
  """ rank values with ties, tie is ranked as the largest rank, now allow nans non-ranked

    Args:
      values(1d np.ma.array): a numeric array

    Returns:
      one vector of asscendant ranks from 1
       ties are kept and label by largest rank 
  """

  V = values
  #print "nonzero(V.mask)=", np.nonzero(V.mask)
  #print "V=", V
  if type(V) == np.ma.MaskedArray:
    #print "V.mask=", V.mask
    nans = (np.nonzero(V.mask)[0]).tolist()      #record indecies
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
  nt = np.nan_to_num(nt)       #filling zeros to nan, shall be no na's from here
  return nt

def noneNormalize(tseries):
  """ no normalizaing

    Args:
      tseries(np.array):  time series matrix

    Returns:
      non normalized tseries
  """

  return np.nan_to_num(tseries)  #filling zeros to nan

def fillMissing(tseries, method): #teseries is 2d matrix unmasked
  """ fill missing data

    Args:
      tseries(np.array):  one time series with no replicates
      method(str):  filling method, choices ['none', 'zero', 'linear', 'slinear', 'nearest', 'quadratic', 'cubic'] 

    Reterns:
      tseries with missing data filled
  """
  
  if method == 'none':
    return tseries                #return with nans
  else:
    y = tseries[np.logical_not(np.isnan(tseries))]
    x = np.array(range(0, len(tseries)))[np.logical_not(np.isnan(tseries))]
    try:
      spline_fit = sp.interpolate.interp1d( x, y, method )
    except:
      print >>sys.stderr, "cannot fill missing values using ", method, "method, fall back to none" 
      return tseries              #return with nans
    yy = np.zeros( len(tseries), dtype='float' )
    for i in range(0, len(tseries)):
      try:
        yy[i] = spline_fit(i)
      except ValueError:
        yy[i] = tseries[i] #keep nans
    return yy
    
def applyAnalysis(firstData, secondData, onDiag=True, delayLimit=3, bootCI=.95, bootNum=1000, pvalueMethod=1000, fTransform=simpleAverage, zNormalize=scoreNormalize):
  """ calculate pairwise LS scores and p-values

    	Args:
    		firstData(np.array): 	numpy data array with correct format factor_num x timespot_num x replicte_num, possibly nans
    		secondData(np.array): 	numpy data array with correct format factor_num x timespot_num x replicte_num, possibly nans
                noDiag(bool):           no results for diagnol comparisions
    		delayLimit(int): 	maximum time unit of delayed response allowed
     		bootCI(float): 		bootstrap confidence interval size, 0 to 1
    		bootNum(int): 		bootstrap number
    		pvalueMethod(int): 	pvalue estimation method and precision
    		ftransform(func): 	summarizing function for replicated data
    		znormalize(func): 	normalizing function for ftransformed data
    		
    	Returns:
    		A LSA table.
    		each row in the table is a list in following format:
    		[ Seq X's Idx, Seq Y's Idx, LS Score, CI_low, CI_high, X's Start Position, 
        	Y's Start Position, Alignment Length, X delay to Y,
        	P-value, Pearson' Correlation, P-value of PCC, Q-value ]
        	
  """	


  firstFactorNum = firstData.shape[0]
  firstRepNum = firstData.shape[1]
  firstSpotNum = firstData.shape[2]
  secondFactorNum = secondData.shape[0]
  secondRepNum = secondData.shape[1]
  secondSpotNum = secondData.shape[2]
  #for now let's assume it is square
  #assert secondFactorNum == firstFactorNum 
  #for now let's assume same rep number
  #assert secondRepNum == firstRepNum 
  #for now let's assume same length
  assert secondSpotNum == firstSpotNum 
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
  timespots = secondSpotNum #same length already assumed
  if pvalueMethod < 0:
    theoPvalue(timespots, delayLimit, 1./np.abs(pvalueMethod), my_decimal)
    #print P_table
  for i in xrange(0, firstFactorNum):
    Xz = np.ma.masked_invalid(firstData[i]) #need to convert to masked array with na's, not F-normalized
    for j in xrange(0, secondFactorNum):
      if onDiag and i<=j:
        continue   #only care lower triangle entries, ignore i=j entries
      Yz = np.ma.masked_invalid(secondData[j])    # need to convert to masked array with na's, not F-normalized
      #if i == 36 or j == 36:
        #print "i=",i,"data=",firstData[i]
        #print "j=",j,"data=",secondData[j]
        #print Xz, Yz
        #print np.all(Yz.mask), np.all(Xz.mask), np.all(Yz.mask) or np.all(Xz.mask)
      if np.all(Yz.mask) or np.all(Xz.mask):  # not any unmasked value in Xz or Yz, all nan in input, continue
        lsaTable[ti] = [i, j, np.nan, np.nan, np.nan, 0, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        pvalues[ti] = np.nan
        pccpvalues[ti] = np.nan
        ti += 1
        continue
      #print "lsa computing..."
      #print "Xz=", Xz
      #print "Yz=", Yz
      LSA_result = singleLSA(Xz, Yz, delayLimit, fTransform, zNormalize, True)                          # do LSA computation
      Smax = LSA_result.score                                                               # get Smax
      Al = len(LSA_result.trace)
      (Xs, Ys, Al) = (LSA_result.trace[Al-1][0], LSA_result.trace[Al-1][1], len(LSA_result.trace))
      #try:
      #except IndexError:
      #  print "Xz=", Xz
      #  print "Yz=", Yz
      #  print "Xs=", Xs, "Ys=", Ys, "Al=", Al  
      #print "bootstrap computing..."
      if bootNum > 0: #do BS
        (Smax, Sl, Su) = bootstrapCI(Xz, Yz, Smax, delayLimit, bootCI, bootNum, fTransform, zNormalize)           # do Bootstrap CI
      else: #skip BS
        (Smax, Sl, Su) = (Smax, Smax, Smax)
      #print "permutation test..."
      #This should be modified.
      if pvalueMethod >= 0:
        lsaP = permuPvalue(Xz, Yz, delayLimit, pvalueMethod, np.abs(Smax), fTransform, zNormalize)          # do Permutation Test
      else:
        #x = np.abs(Smax)*np.sqrt(timespots) # x=Rn/sqrt(n)=Smax*sqrt(n)
        lsaP = readPvalue(np.abs(Smax), timespots, my_decimal) #one-tailed  
      pvalues[ti] = lsaP
      #print "PPC..." 
      #print np.mean(Xz, axis=0), np.mean(Yz, axis=0)
      (PCC, P_PCC) = sp.stats.pearsonr(np.ma.average(Xz, axis=0), np.ma.average(Yz, axis=0)) # two tailed p-value
      #P_PCC = P_PCC/2   # one tailed p-value
      #(DPCC, P_DPCC) = sp.stats.pearsonr(np.ma.average(Xz[:,Xs-1:Xs+Al],
      #print Xs, Ys, Al, Ys-Xs
      if Xs <= Ys:
        #print Xz[:,:Al].shape
        #print Yz[:,(Ys-Xs):(Ys-Xs)+Al].shape
        (SPCC, P_SPCC) = sp.stats.pearsonr(np.ma.average(Xz[:,:Al], axis=0), np.ma.average(Yz[:,(Ys-Xs):(Ys-Xs+Al)], axis=0)) # corr for shifted-cut seq
      else:
        #print Xz[:,(Xs-Ys):(Xs-Ys)+Al].shape
        #print Yz[:,:Al].shape
        (SPCC, P_SPCC) = sp.stats.pearsonr(np.ma.average(Xz[:,(Xs-Ys):(Xs-Ys+Al)], axis=0), np.ma.average(Yz[:,:Al], axis=0)) # corr for shifted-cut seq

      # need +epsilon to avoid all zeros
      pccpvalues[ti] = P_PCC
      lsaTable[ti] = [i, j, Smax, Sl, Su, Xs, Ys, Al, Xs-Ys, lsaP, PCC, P_PCC, SPCC, P_SPCC]
      ti += 1
      #print "finalizing..."

  #print lsaTable
  #print "qvalues ..."
  #pvalues = np.ma.masked_invalid(pvalues)
  #pccpvalues = np.ma.masked_invalid(pccpvalues)
  #print "pvalues", pvalues, np.isnan(np.sum(pvalues))
  #print "pccpvalues", pccpvalues, np.isnan(np.sum(pccpvalues))
  qvalues = storeyQvalue( pvalues )
  pccqvalues = storeyQvalue( pccpvalues )
  #print "qvalues", qvalues
  #print "pccqvalues", pccqvalues

  for i in xrange(0, len(qvalues)):
    lsaTable[i].append( qvalues[i] )
    lsaTable[i].append( pccqvalues[i] )

  #print lsaTable
  return lsaTable

def test():
  """ self test script
  """
  np.seterr(all='raise')
  print >>sys.stderr, "###testing###"
  test_data = np.array( [ [ [ 3, 2, 3, 4], [1, 6, 2, 3], [6, 4, 6, 8], [np.nan, np.nan, np.nan, np.nan] ], 
                          [ [np.nan, 2, np.nan, 3], [4, 5, 3, np.nan], [1, 1, 1, 1], [ np.nan, np.nan, np.nan, np.nan ] ] ], dtype='float' )
  test_fN = test_data.shape[0]
  test_rN = test_data.shape[1]
  test_tN = test_data.shape[2]
  print >>sys.stderr, "fN, tN, rN", test_fN, test_tN, test_rN
  print >>sys.stderr, "test_data", test_data
  print >>sys.stderr, "---fillMissing---"
  print >>sys.stderr, "none:", test_data[1][0], "->", fillMissing(test_data[1][0], 'none')
  print >>sys.stderr, "zero (spline at zero order):", test_data[1][0], "->", fillMissing(test_data[1][0], 'zero')
  print >>sys.stderr, "linear:", test_data[1][0], "->", fillMissing(test_data[1][0], 'linear')
  print >>sys.stderr, "quadratic:", test_data[1][0], "->", fillMissing(test_data[1][0], 'quadratic')
  print >>sys.stderr, "cubic:", test_data[1][0], "->", fillMissing(test_data[1][0], 'cubic')
  print >>sys.stderr, "slinear:", test_data[1][0], "->", fillMissing(test_data[1][0], 'slinear')
  print >>sys.stderr, "nearest:", test_data[1][0], "->", fillMissing(test_data[1][0], 'nearest')
  print >>sys.stderr, "###use none as clean data###"
  clean_data = test_data
  for i in xrange(0, test_fN):
    for j in xrange(0, test_rN):
      clean_data[i][j] = fillMissing(test_data[i][j], 'none')
  print >>sys.stderr, "clean_data", clean_data
  print >>sys.stderr, "---masked_array---"
  masked_data = np.ma.masked_invalid(clean_data)
  print >>sys.stderr, "masked_data", masked_data
  print >>sys.stderr, "---tied-rank---"
  print >>sys.stderr, "data:", masked_data[1][0], "->", tied_rank(masked_data[1][0])
  print >>sys.stderr, "data:", masked_data[1][1], "->", tied_rank(masked_data[1][1])
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
  pvalues = np.array([0.01, 0.2, 0.03, 0.4, 0.05, np.nan, 0.03, 0.4, 0.03, 0.3], dtype='float')
  print >>sys.stderr, "pvalues:", pvalues 
  print >>sys.stderr, "qvalues:", storeyQvalue(pvalues)
  pvalues = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.02, 0.03, 0.04, 0.03, 0.03], dtype='float')
  print >>sys.stderr, "pvalues:", pvalues 
  print >>sys.stderr, "qvalues:", storeyQvalue(pvalues)
  print >>sys.stderr, "---singleLSA---"
  print >>sys.stderr, "input data:", scoreNormalize(simpleAverage(masked_data[0])), \
                                     scoreNormalize(simpleAverage(masked_data[1]))
  lsar=singleLSA(masked_data[0], masked_data[1], delayLimit=1, fTransform=simpleAverage, zNormalize=scoreNormalize, keepTrace=True)
  print >>sys.stderr, "lsar.score=", lsar.score 
  Al=len(lsar.trace)
  print >>sys.stderr, "lsar.align=",(lsar.trace[Al-1][0], lsar.trace[Al-1][1], len(lsar.trace)) 
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
  print >>sys.stderr, applyAnalysis(clean_data, clean_data, True, 1, .95, 10, 10, sdAverage, scoreNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, clean_data, True, 1, .95, 10, 10, simpleAverage, scoreNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, clean_data, True, 1, .95, 10, 10, simpleMedian, scoreNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, clean_data, True, 1, .95, 10, 10, madMedian, scoreNormalize)

if __name__=="__main__":
  print "hello world!"
  test()
  exit(0)


#############  FUNCTIONS NO IN USE AFTER THIS LINE ####################

