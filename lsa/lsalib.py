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

#Author Note:
#Clearification of current Handling of Data Matrix:
# -input-> numpy.ndarray(with na) -fillMissing-> np.ndarray(with na) -applyAnalysis->
# np.masked(na masked) -fTransform-> np.masked(na masked) -zNormalize-> np.ndarray(masked=0, no na)

#Considering not using numpy.masked array (their implementation having errors, how to report?)
#Considering using R for simple numerics, rpy or use swig+R?

#import public resources
import csv, sys, os, random
import numpy as np
#import numpy.ma as np.ma
import scipy as sp
import scipy.interpolate
import scipy.stats
#R through Rpy
import rpy2.rlike.container as rlc
import rpy2.robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
ro.conversion.py2ri = numpy2ri
r = ro.r

#print '''setwd("%s")''' % os.environ.get('PWD')
r('''setwd("%s")''' % os.environ.get('PWD'))
r('''options(warn=-1)''') 

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
    #Xz = percentileNormalize( X )												#N: normalize function
    #S_max = singleLSA( X, Y, D )											#D: delayLimit	
    #(S_low, S_up) = bootstrapCI( X, Y, D, alpha, B )						#B: bootNum, alpha: interval size
    #p = permuPvalue( X, Y, S_max, P )										#P: pvalueMethod
    #q = storeyQvalue( p_set )												#p_set: set of p-values
    #lsaTable = applyAnalysis( rawData, delayLimit=D, ftransform=F, znormalize=N, bootCI=alpha, bootNum=B, pvalueMethod=P )  
"""

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

###############################
# applyAnalsys
# return a list-of-list, where the i-th list is:
#	[ f1, f2, L-score, L_low, L_up, X_start, Y_start, length, P/N, P-value, PCC,  PCC P-val,  Qvalue  ]
#	[ 1,  2,  3,       4,     5,    6,       7,       8,      9,   10,      11,   12,         13      ]
############################### 

def rpy_spearmanr(Xz, Yz):
  sr=r('''cor.test''')(Xz,Yz,method='spearman')
  return (sr[3][0],sr[2][0])

def calc_spearmanr(Xz, Yz, sfunc=rpy_spearmanr):
  mask = np.logical_or(Xz.mask, Yz.mask)
  Xz.mask = mask
  Yz.mask = mask
  (SCC, P_SCC) = sfunc(Xz.compressed(), Yz.compressed()) # two tailed p-value
  return (SCC, P_SCC)

def calc_pearsonr(Xz, Yz, pfunc=scipy.stats.pearsonr):
  mask = np.logical_or(Xz.mask, Yz.mask)
  Xz.mask = mask
  Yz.mask = mask
  (PCC, P_PCC) = pfunc(Xz.compressed(), Yz.compressed()) # two tailed p-value
  return (PCC, P_PCC)

def calc_shift_corr(Xz, Yz, D, corfunc=calc_pearsonr): 
  d_max=0
  r_max=0
  p_max=1
  #print "Xz=", Xz
  #print "Yz=", Yz
  for d in range(-D, D+1):
    # i = Xs-Ys
    #print "d=", d
    if d < 0:
      X_seg = Xz[:(len(Xz)+d)]
      Y_seg = Yz[-d:len(Yz)]
    elif d == 0:
      X_seg = Xz
      Y_seg = Yz
    else:
      X_seg = Xz[d:len(Xz)]
      Y_seg = Yz[:len(Yz)-d]
    #print "Xz=", X_seg.shape
    #print "Yz=", Y_seg.shape
    assert len(X_seg) == len(Y_seg)
    mask = np.logical_or(X_seg.mask, Y_seg.mask)
    X_seg.mask = mask
    Y_seg.mask = mask
    cor = corfunc(X_seg, Y_seg)
    if np.abs(cor[0]) >= np.abs(r_max):
      r_max = cor[0] 
      d_max = d
      p_max = cor[1]
  return (r_max, p_max, d_max)
        
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
  
  #print "f1=", isinstance(series1,np.ma.core.MaskedArray)
  #print "f2=", isinstance(series2,np.ma.core.MaskedArray)
  #print "f1=", isinstance(fTransform(series1),np.ma.core.MaskedArray)
  #print "f2=", isinstance(fTransform(series2),np.ma.core.MaskedArray)
  #try:
  lsad=compcore.LSA_Data(delayLimit, zNormalize(fTransform(series1)), zNormalize(fTransform(series2)))
  #except NotImplementedError:
  #  print series1, series1.mask, series2, series2.mask
  #  print fTransform, fTransform(series1), fTransform(series1).mask, fTransform(series2), fTransform(series2).mask
  #  print zNormalize, zNormalize(fTransform(series1)), zNormalize(fTransform(series1)).mask, zNormalize(fTransform(series2)), zNormalize(fTransform(series2)).mask
  #  quit()
  lsar=compcore.DP_lsa(lsad, keepTrace)
  del lsad
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

def ma_average(ts, axis=0):
  #print "before ma_average", ts.shape
  ns = np.ma.mean(ts, axis=0)
  #print "after ma_average", ns.shape
  if type(ns.mask) == np.bool_:       #fix broken ma.mean, mask=False instead of [False, ...] for all mask
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

def readPvalue(P_table, R, N, x_sd=1, M=1, alpha=1, beta=1, x_decimal=3):  
  # R=observed range, N=timepoints, x_sd=std.dev of single series, M=replicates, alpha=1-portion of zero in X, beta=1-portion of zero in Y
  # x' = R*M/(alpha*beta*sqrt(N)*sd) * 10^(x_decimal)
  # has to ceil the x value to avoid round to 0, which is not amenable to calculation
  try:
    xi = int(np.around(R*M/np.sqrt(alpha*beta*N)*(10**x_decimal)))    #
  except OverflowError:
    #print alpha, beta
    #print "x=", np.around(R*M/np.sqrt(alpha*beta*N)*(10**x_decimal))
    #quit()
    return np.nan
  if xi in P_table:
    return P_table[xi]
  else:
    return 0.

def theoPvalue(Rmax, Dmax=0, precision=.001, x_decimal=2):   #let's produce 2 tail-ed p-value
  # the produced P_table is for P(X>=x) when X=(R(D)/sqrt(n)) and indexed by xi=x*(10^x_decimal) 
  # so, x=xi/(10^x_decimal)
  # no read for recaculating P_table when read proper scaled x' = R*M/(alpha*beta*sqrt(N)*sd) * 10^(x_decimal)
  # print np.linspace(0,timespots,timespots/10**(-x_decimal)+1) #x_decimal is for augment x_index
  Rmax = np.max((Rmax, Rmax_min))
  P_table = dict()
  for xi in xrange(0,Rmax*10**(x_decimal)+1): 
    if xi == 0:
      P_table[xi] = 1
      continue
    #print x
    x = xi/float(10**(x_decimal)) #standard x with variance corrected
    xx = x**2
    pipi_over_xx = pipi/xx
    alpha = precision
    B = 2*Dmax+1
    #print pipi_over_xx
    #print alpha
    #print np.log(alpha*xx*(1-np.exp(-pipi_over_xx))/(8**B)/2)
    #print np.log(alpha*xx*(1-np.exp(-pipi_over_xx))/(8**B)/2) / pipi_over_xx
    Kcut = np.max((kcut_min, int(np.ceil( .5 - np.log( (alpha/(2**B-1))**(1/B) *xx*(1-np.exp(-pipi_over_xx))/8/2 )/pipi_over_xx ))))
    #Kcut = 200
    A = 1/xx
    Rcdf = 0     # root of cdf

    for k in xrange(1,Kcut+1):
      C = (2*k-1)**2
      Rcdf = Rcdf + (A+pipi_inv/C)*np.exp(-C*pipi_over_xx/2)
      P_two_tail = 1 - (8**B)*(Rcdf**B)
      #if x==2.52:
      #  print "k=",k, "A=",A, "B=",B, "C=",C, "F=",(8**B)*(R**B), "dR=",(A+pipi_inv/C)*np.exp(-C*pipi_over_xx/2), "P=",P,"pipi_inv=",pipi_inv,"pipi_over_xx=",pipi_over_xx 
    #print xi, x, Kcut, P, P/2;

    #print "xi=", xi, "Kut=", Kcut, "k=", k
    #if P_two_tail <= precision and Kcut != :
    #  P_table[xi] = precision
    #  continue
    #else:
    P_table[xi] = P_two_tail #return two tailed probability

  return P_table
	
def permuPvalue(series1, series2, delayLimit, precisionP, Smax, fTransform, zNormalize):
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
  PP_set = np.zeros(precisionP, dtype='float')
  Xz = zNormalize(fTransform(series1))
  Y = np.ma.array(series2)                                               #use = only assigns reference, must use a constructor
  for i in xrange(0, precisionP):
    np.random.shuffle(Y.T)
    lsad.assign( delayLimit, Xz, zNormalize(fTransform(Y)) )
    PP_set[i] = compcore.DP_lsa(lsad, False).score
  #PP_set[pvalueMethod]=Smax                                               #the original test shall not be considerred
  #print "PP_set", PP_set, PP_set >= Smax, np.sum(PP_set>=Smax), np.float(pvalueMethod)
  if Smax >= 0:
    P_two_tail = np.sum(np.abs(PP_set) >= Smax)/np.float(precisionP)
  else:
    P_two_tail = np.sum(-np.abs(PP_set) <= Smax)/np.float(precisionP)
  return P_two_tail

#Q_lam_step = 0.05
#Q_lam_max = 0.95
#qvalue(p=NULL, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=NULL, robust=FALSE, gui=FALSE, 
#       smooth.df=3, smooth.log.pi0=FALSE)
def R_Qvalue(pvalues, lam=np.arange(0,Q_lam_max,Q_lam_step), method='smoother', robust=False, smooth_df=3):
  try:
    pvalues_not_nan = np.logical_not(np.isnan(pvalues))
    pvalues_input = pvalues[pvalues_not_nan] 
    r('''library(qvalue)''')
    #print "pvalues=", pvalues
    #print "pvalues_input=", pvalues_input
    qvalues=r['''qvalue'''](p=pvalues_input, **{'lambda':lam, 'pi0.method':method, 'robust':robust, 'smooth.df':3})
    qvalues_return = [np.nan]*len(pvalues)
    #print "qvalues=", qvalues
    #print len(qvalues)
    #print qvalues[2] -- this is $qvalues
    j = 0
    for i in xrange(0, len(pvalues)):
      if not np.isnan(pvalues[i]):
        qvalues_return[i]=qvalues[2][j]   #second item is calculated qvalues
        j = j+1
  except:
    #print >>sys.stderr, "caution: q-value estimation error"
    print >>sys.stderr, "from R: unusable pvalues -> ", pvalues_input
    qvalues_return=[np.nan]*len(pvalues)
  #print "qvalues_return=", qvalues_return
  return qvalues_return

def storeyQvalue(pvalues, lam=np.arange(0,Q_lam_max,Q_lam_step), method='smoother', robust=False, smooth_df=3):
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
  
  try:
    #assume pvalues is an array with possible nans
    #get a list of non nan pvalues, do the same procedure, putback to original list
    mpvalues = np.ma.masked_invalid(pvalues,copy=True)
    rpvalues = mpvalues[~mpvalues.mask]
    p_num = len(pvalues)
    rp_num = len(rpvalues)
    #rp_nnz = len(np.nonzero(rpvalues))
    #print "p_num=", p_num, "rp_num=", rp_num

    if rp_num <= 1:
      #print >>sys.stderr, "WARN: not enough number of pvalues for q-value evaluation! nan will be filled!"
      return np.array( [np.nan] * p_num, dtype='float')

    rp_max = np.max(rpvalues)
    rp_lam = lam[lam<rp_max]

    if len(rp_lam) <= 1:
      return np.array( [np.nan if np.isnan(pvalues[i]) else 0 for i in xrange(0,p_num)],  dtype='float')

    #print "rpvalues=", rpvalues, rp_num, rp_lam

    pi_set = np.zeros(len(rp_lam), dtype='float')
    for i in xrange(0, len(rp_lam)): 
      pi_set[i] = np.mean(rpvalues>=rp_lam[i])/(1-rp_lam[i])

    #print "p=",pvalues
    #print len(pi_set), rp_max, len(rp_lam) #rp_nnz

    #print "pi_set=", pi_set
    if method=='smoother':
      spline_fit = sp.interpolate.interp1d(rp_lam, pi_set, kind=smooth_df)
      pi_0 = spline_fit(np.max(rp_lam))
      #print "pi_0=", pi_0
      pi_0 = np.max( [np.min( [np.min(pi_0), 1]), 0] ) #0<=pi_0<=1
      #print "pi_0=", pi_0
      if pi_0 == 0:
        #print >>sys.stderr, "WARN: smoother method not working, fall back to bootstrap"
        #print pi_set, rp_max, rp_lam, pi_0
        method='bootstrap'
    #print "pi_set=", pi_set
    #print "pi_0=", pi_0

    if method=='bootstrap':                            #bootstrap
      pi_min = np.min(pi_set)
      mse = np.zeros((100, len(rp_lam)), dtype='float')
      pi_set_boot = np.zeros((100, len(rp_lam)), dtype='float')
      for j in xrange(0, 100):
        p_boot = sample_wr(rpvalues, rp_num)
        for i in xrange(0, len(rp_lam)):
          pi_set_boot[j][i] = np.mean(p_boot>=rp_lam[i])/(1-rp_lam[i]) 
        mse[j] = (pi_set_boot[j]-pi_min)**2
      min_mse_j = np.argmin(mse)
      pi_0 = np.min(pi_set_boot[min_mse_j])
      #print "pi_0=", pi_0
      pi_0 = np.max([np.min( [np.min(pi_0), 1]), ]) #0<=pi_0<=1
      #print "pi_0=", pi_0
      if pi_0 == 0:
        #print >>sys.stderr, "WARN: bootstrap method not working, cannot estimate qvalues"
        pi_0 = Q_lam_step #non-fixable peculiar case, show some reasonable results, refer to Tibshirani et al.
        #return np.array( [np.nan] * p_num, dtype='float' )

    #print "pi_0=", pi_0
    rp_argsort = np.argsort(rpvalues)                     #np.nan will be sorted as maximum (the largest rank)
    #print "argsort of rps=", rpvalues[rp_argsort]
    rp_ranks = tied_rank(rpvalues)
    #print "tied rank of rps=", rp_ranks
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
 
    #if np.all(np.isnan(qvalues)):
    #print method
    #print "q=",qvalues
    #print "rp_sort",rp_argsort
    #print "rp_rank",rp_ranks
    #print "qs=", qvalues
  except:
    print >>sys.stderr, "caution: q-value estimation error"
    print >>sys.stderr, "from scipy: unusable pvalues -> ", rpvalues
    qvalues=np.array( [np.nan] * p_num, dtype='float')

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
  Xf = ma_average(tseries, axis=0)
  return Xf

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
    return simpleAverage(tseries)                       #sd = 0, fall back to simpleAverage
  if np.any(sd.mask) or (np.ma.sum(sd==0))>0:
    #print sd, sd.mask, sd==0
    return simpleAverage(tseries)                       #sd = 0, fall back to simpleAverage
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
  Xf = ma_average(tseries, axis=0)*(1/sd)*(1/np.ma.sum(1/sd))*(1/sd)   #sd-weighted sample
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
    return simpleMedian(tseries)                  #mad = 0, fall back to simpleMedian
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
  assert type(values) == np.ma.MaskedArray
  #print "values=", values
  #print "values.mask=", values.mask
  V = np.ma.asarray(values)
    #print "nonzero(V.mask)=", np.nonzero(V.mask)
    #print "V=", V
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
  sV = np.array( [ (2*v_cum[V[i]]-v_num[V[i]]+1)/2 for i in xrange(0, len(V)) ], dtype='float' ) #sorted V, break tie by average
  #print "sV=", sV
  #print "nans=", nans
  
  #if type(values) == np.ma.MaskedArray:
    #mV = sV
  for idx in nans:
    sV = np.insert(sV, idx, np.nan)   #insert nan to original place, need return masked? yes
  sV = np.ma.masked_invalid(sV,copy=True)
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

def percentileNormalize(tseries):
  """ score normalizing

    Args:
      tseries(np.array): 1-d time series
    
    Returns:
      score normalized time series
  """
  ranks = tied_rank(tseries)
  #print "ranks=", ranks
  nt = np.ma.masked_invalid(sp.stats.distributions.norm.ppf( ranks/(len(ranks)-np.sum(ranks.mask)+1) ),copy=True)
  #print "nt=", nt
  #nt = np.nan_to_num(nt)              #filling zeros to nan, shall be no na's from here on
  nt = nt.filled(fill_value=0)         #filling zeros to nan, shall be no na's from here on
  return nt

def percentileZNormalize(tseries):
  """ score normalizing

    Args:
      tseries(np.array): 1-d time series
    
    Returns:
      score normalized time series
  """
  ranks = tied_rank(tseries)
  #print "tseries=", tseries
  #print "ranks=", ranks
  #print "ranks.mask", ranks.mask
  nt = np.ma.masked_invalid(sp.stats.distributions.norm.ppf( ranks/(len(ranks)-np.sum(ranks.mask)+1) ),copy=True)
  try:
    zt = (nt - np.ma.mean(nt, axis=0))/np.ma.std(nt)
  except FloatingPointError:
    zt = nt - np.ma.mean(nt, axis=0)
  #print "nt=", nt
  #nt = np.nan_to_num(nt)              #filling zeros to nan, shall be no na's from here on
  zt = zt.filled(fill_value=0)         #filling zeros to nan, shall be no na's from here on
  return zt

def noZeroNormalize(tseries):
  """ score normalizing

    Args:
      tseries(np.array): 1-d time series
    
    Returns:
      score normalized time series
  """
  #print "with zero tseries=", tseries, tseries.mask
  #print "t=", tseries
  nt = np.ma.masked_equal(tseries, 0)
  #print "nt=", nt, nt.mask
  if type(nt.mask) == np.bool_:       #fix broken ma.mean, mask=False instead of [False, ...] for all mask
    nt.mask = [nt.mask] * nt.shape[0]
  #print "none zero tseries=", nt, nt.mask
  ranks = tied_rank(nt)
  #print ranks, ranks.mask, len(ranks)-np.sum(ranks.mask), ranks/(np.sum(ranks.mask)+1) 
  nt = np.ma.masked_invalid(sp.stats.distributions.norm.ppf( ranks/(len(ranks)-np.sum(ranks.mask)+1) ),copy=True)
  #print np.ma.mean(nt, axis=0), np.ma.std(nt, axis=0)
  #print (nt - np.ma.mean(nt, axis=0))*(1/np.ma.std(nt, axis=0))
  try:
    zt = (nt - np.ma.mean(nt, axis=0))*(1/np.ma.std(nt, axis=0))
  except FloatingPointError:
    zt = nt - np.ma.mean(nt, axis=0)
  #nt = np.nan_to_num(nt)              #filling zeros to nan, shall be no na's from here on
  #print "nt=", nt
  zt = zt.filled(fill_value=0)         #filling zeros to nan, shall be no na's from here on
  #print "zt=", zt
  return zt

def noneNormalize(tseries):
  """ no normalizaing

    Args:
      tseries(np.array):  time series matrix

    Returns:
      non normalized tseries
  """

  nt = tseries.filled(fill_value=0)   #filling zeros to nan, shall be no na's from here on
  return nt

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
    x = np.array(range(0, len(tseries)), dtype='float')[np.logical_not(np.isnan(tseries))]
    try:
      spline_fit = sp.interpolate.interp1d( x, y, method )
    except:
      #print >>sys.stderr, "cannot fill missing values using ", method, "method, fall back to none" 
      return tseries              #return with nans
    yy = np.zeros( len(tseries), dtype='float' )
    for i in range(0, len(tseries)):
      if not np.isnan(tseries[i]):
        yy[i] = tseries[i]
      else:
        try:
          yy[i] = spline_fit(i)
        except ValueError:
          yy[i] = tseries[i] #keep nans

    return yy
    
def applyAnalysis(firstData, secondData, onDiag=True, delayLimit=3, minOccur=.5, bootCI=.95, bootNum=0, pvalueMethod='perm', precisionP=1000,\
    fTransform=simpleAverage, zNormalize=noZeroNormalize, varianceX=1, resultFile=open("tmp.lsa","w"), \
    firstFactorLabels=None, secondFactorLabels=None, qvalueMethod='R'):
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

  col_labels= ['X','Y','LS','lowCI','upCI','Xs','Ys','Len','Delay','P','PCC','Ppcc','SPCC','Pspcc','Dspcc','SCC','Pscc','SSCC','Psscc','Dsscc',
            'Q','Qpcc','Qspcc','Qscc','Qsscc','Xi','Yi']
  print >>resultFile,  "\t".join(col_labels)

  firstFactorNum = firstData.shape[0]
  firstRepNum = firstData.shape[1]
  firstSpotNum = firstData.shape[2]
  secondFactorNum = secondData.shape[0]
  secondRepNum = secondData.shape[1]
  secondSpotNum = secondData.shape[2]
  if not firstFactorLabels:
    firstFactorLabels= [str(v) for v in range(1, firstFactorNum+1)]
  if not secondFactorLabels:
    secondFactorLabels= [str(v) for v in range(1, secondFactorNum+1)]
  #for now let's assume same rep number
  #for now let's assume same length
  assert secondSpotNum == firstSpotNum 
  assert secondRepNum == firstRepNum
  if onDiag:  # if assigned jobs are on the diagnal
    assert firstFactorNum == secondFactorNum
    pairwiseNum = firstFactorNum*(firstFactorNum-1)/2
  else:
    pairwiseNum = firstFactorNum*secondFactorNum
  lsaTable = [None]*pairwiseNum
  pvalues = np.zeros(pairwiseNum, dtype='float')
  pccpvalues = np.zeros(pairwiseNum, dtype='float')
  spccpvalues = np.zeros(pairwiseNum, dtype='float')
  sccpvalues = np.zeros(pairwiseNum, dtype='float')
  ssccpvalues = np.zeros(pairwiseNum, dtype='float')
  #print factorNum, repNum, spotNum, lsaTable, pvalues

  timespots = secondSpotNum #same length already assumed
  replicates = firstRepNum
  stdX = np.sqrt(varianceX) #make comparable with previous command line
  ti = 0
  
  if qvalueMethod in ['scipy']:
    qvalue_func = storeyQvalue
  else:
    qvalue_func = R_Qvalue 

  if pvalueMethod in ['theo','mix']:
    #P_table = theoPvalue(D=0, precision=.0001, x_decimal=3)   #let's produce 2 tail-ed p-value
    P_table = theoPvalue(Rmax=timespots, Dmax=delayLimit, precision=1./precisionP, x_decimal=my_decimal)
    #print P_table
  for i in xrange(0, firstFactorNum):
    Xz = np.ma.masked_invalid(firstData[i], copy=True) #need to convert to masked array with na's, not F-normalized
    Xz_badOccur = np.sum(np.logical_not(np.isnan(ma_average(Xz)), ma_average(Xz)==0))/float(timespots) < minOccur
    if Xz.shape[1] == None: #For 1-d array, convert to 2-d
      Xz.shape = (1, Xz.shape[0])
    for j in xrange(0, secondFactorNum):
      if onDiag and i>=j:
        continue   #only care lower triangle entries, ignore i=j entries
      Yz = np.ma.masked_invalid(secondData[j], copy=True)    # need to convert to masked array with na's, not F-normalized
      Yz_badOccur = np.sum(np.logical_not(np.isnan(ma_average(Yz)), ma_average(Yz)==0))/float(timespots) < minOccur
      #print i+1, np.sum(np.logical_not(np.isnan(ma_average(Xz)), ma_average(Xz)==0))/float(timespots), Xz_badOccur, \
      #	    j+1, np.sum(np.logical_not(np.isnan(ma_average(Yz)), ma_average(Yz)==0))/float(timespots), Yz_badOccur
      #if Yz_badOccur:
	#print ma_average(Yz), ma_average(Yz).mask, quit()
      if Yz.shape[1] == None: #For 1-d array, convert to 2-d
        Yz.shape = (1, Yz.shape[0])
      #if i == 36 or j == 36:
        #print "i=",i,"data=",firstData[i]
        #print "j=",j,"data=",secondData[j]
        #print Xz, Yz
        #print np.all(Yz.mask), np.all(Xz.mask), np.all(Yz.mask) or np.all(Xz.mask)

      #control minZeroPercent or allNA here
      
      if np.all(Yz.mask) or np.all(Xz.mask) or Xz_badOccur or Yz_badOccur:  # not any unmasked value in Xz or Yz, all nan in input, warn code -1
      #lsaTable[ti] = [i, j, Smax,   Sl,     Su,     Xs,Ys,Al,Xs-Ys, lsaP,   PCC,    P_PCC,  SPCC,   P_SPCC, D_SPCC, SCC,    P_SCC,  SSCC,   P_SSCC, D_SSCC]
        lsaTable[ti] =[i, j, np.nan, np.nan, np.nan, -1, -1, -1, -1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        pvalues[ti] = np.nan
        pccpvalues[ti] = np.nan
        spccpvalues[ti] = np.nan
        sccpvalues[ti] = np.nan
        ssccpvalues[ti] = np.nan
        ti += 1
        continue
      #print "lsa computing..."
      #try:
      #print "Xf%d="%i, Xz
      #print "Xz%d="%i, zNormalize(fTransform(Xz))
      #print "Xf%d="%j, Yz
      #print "Xz%d="%j, zNormalize(fTransform(Yz))
      LSA_result = singleLSA(Xz, Yz, delayLimit, fTransform, zNormalize, True)                          # do LSA computation
      #except NotImplementedError:
      #  print "Xz=", Xz
      #  print "Yz=", Yz
      #  quit()
          
      Smax = LSA_result.score                                                               # get Smax
      Al = len(LSA_result.trace)
      (PCC, P_PCC) = calc_pearsonr(ma_average(Xz, axis=0), ma_average(Yz, axis=0)) # two tailed p-value
      (SCC, P_SCC) = calc_spearmanr(ma_average(Xz, axis=0), ma_average(Yz, axis=0))
      
      try:
        (SPCC, P_SPCC, D_SPCC) = calc_shift_corr(ma_average(Xz, axis=0), ma_average(Yz, axis=0), delayLimit, calc_pearsonr) 
        # corr for shifted-cut seq
        (SSCC, P_SSCC, D_SSCC) = calc_shift_corr(ma_average(Xz, axis=0), ma_average(Yz, axis=0), delayLimit, calc_spearmanr) 
        # corr for shifted-cut seq
        #if np.isnan(SSCC) or np.isnan(SPCC): 
        #  print "Al=", Al, "X shape", Xz.shape, "Y shape", Yz.shape
        #  print "Xs=", Xs, X_seg.shape
        #  print "Ys=", Ys, Y_seg.shape
        #  quit()
      except FloatingPointError:
        (SPCC, P_SPCC, D_SPCC) = (np.nan, np.nan, np.nan)
        (SSCC, P_SSCC, D_SSCC) = (np.nan, np.nan, np.nan)
        #print np.ma.mean(Xz[:,:Al], axis=0), np.ma.mean(Xz[:,:Al], axis=0).mask, \
        #  np.ma.mean(Yz[:,(Ys-Xs):(Ys-Xs+Al)], axis=0), np.ma.mean(Yz[:,(Ys-Xs):(Ys-Xs+Al)], axis=0).mask
        #quit()
      pccpvalues[ti] = P_PCC
      spccpvalues[ti] = P_SPCC
      try:
        sccpvalues[ti] = P_SCC
      except ValueError:
        #print "error at lsalib P_SCC", P_SCC  #sometimes produce [], why?
        sccpvalues[ti] = np.nan
      ssccpvalues[ti] = P_SSCC

      if Al == 0: #handel align impossibility, usually too many nas'
        #(SCC, P_SCC) = calc_spearmanr(ma_average(Xz, axis=0), ma_average(Yz, axis=0), statlib.stats.lspearmanr)
        pvalues[ti] = np.nan
        #lsaTable[ti] = [i, j, Smax,   Sl,     Su,     Xs,Ys,Al,Xs-Ys, lsaP,   PCC, P_PCC,  SPCC,   P_SPCC, D_SPCC, SCC, P_SCC,  SSCC, P_SSCC, D_SSCC]
        lsaTable[ti] =  [i, j, np.nan, np.nan, np.nan, 0, 0, 0, 0,     np.nan, PCC, P_PCC,  SPCC, P_SPCC, D_SPCC, SCC, P_SCC, SSCC, P_SSCC, D_SSCC]
        ti += 1
        continue

      (Xs, Ys, Al) = (LSA_result.trace[Al-1][0], LSA_result.trace[Al-1][1], len(LSA_result.trace))
      #try:
      #  (Xs, Ys, Al) = (LSA_result.trace[Al-1][0], LSA_result.trace[Al-1][1], len(LSA_result.trace))
      #except IndexError:
      #  print "Xz=", Xz
      #  print "Yz=", Yz
      #  print "Xs=", Xs, "Ys=", Ys, "Al=", Al  
      #  quit()
      #print "PPC..." 
      #print np.mean(Xz, axis=0), np.mean(Yz, axis=0)

      #if Xs <= Ys:
        #print Xz[:,:Al].shape
        #print Yz[:,(Ys-Xs):(Ys-Xs)+Al].shape
      #  X_seg = Xz[:,:(Xz.shape[1]-(Ys-Xs))]
      #  Y_seg = Yz[:,(Ys-Xs):Yz.shape[1]]
      #else:
      #  X_seg = Xz[:,(Xs-Ys):Xz.shape[1]]
      #  Y_seg = Yz[:,:(Yz.shape[1]-(Xs-Ys))]
      #print "Al=", Al, "X shape", Xz.shape, "Y shape", Yz.shape
      #print "Xs=", Xs, X_seg.shape
      #print "Ys=", Ys, Y_seg.shape
      #if Xs != Ys:
      #  quit()


      #This Part Must Follow Static Calculation Part
      #Otherwise Yz may be changed, now it is copied
      #np.ma.array(copy=True) to copy, otherwise is only reference
      promisingP = 0.05 #this is according to convention of P-value
      if pvalueMethod in ['theo', 'mix']:
        #x = np.abs(Smax)*np.sqrt(timespots) # x=Rn/sqrt(n)=Smax*sqrt(n)
        #alpha=1-{#(nan+zeros)/timespots}
        Xz_alpha = 1 - np.sum(zNormalize(fTransform(Xz))==0)/float(timespots)
        Yz_beta = 1 - np.sum(zNormalize(fTransform(Yz))==0)/float(timespots)
        if Xz_alpha==0 or Yz_beta==0:
          lsaP = np.nan
        else:
          lsaP = readPvalue(P_table, R=np.abs(Smax)*timespots, N=timespots, x_sd=stdX, M=replicates, alpha=Xz_alpha, beta=Yz_beta, x_decimal=my_decimal)

      if (pvalueMethod in ['mix'] and lsaP<=promisingP) or (pvalueMethod in ['perm']):
        Xp = np.ma.array(Xz,copy=True)
        Yp = np.ma.array(Yz,copy=True)
        lsaP = permuPvalue(Xp, Yp, delayLimit, precisionP, np.abs(Smax), fTransform, zNormalize)          # do Permutation Test

      pvalues[ti] = lsaP
      #print "bootstrap computing..."
      if bootNum > 0: #do BS
        Xb = np.ma.array(Xz,copy=True)
        Yb = np.ma.array(Yz,copy=True)
        (Smax, Sl, Su) = bootstrapCI(Xb, Yb, Smax, delayLimit, bootCI, bootNum, fTransform, zNormalize)           # do Bootstrap CI
      else: #skip BS
        (Smax, Sl, Su) = (Smax, Smax, Smax)

      #uncomment to verify Xz, Yz unchanged by bootstrap and permutation
      #print 'i=',i,"Data=",firstData[i],"Xz=",Xz,"mask=",Xz.mask, "mask_invalid=", np.ma.masked_invalid(firstData[i])
      #print 'j=',j,"Data=",secondData[j],"Yz=",Yz,"mask=",Yz.mask, "mask_invalid=", np.ma.masked_invalid(secondData[j])
      #print ma_average(Xz, axis=0), "mask=", ma_average(Xz, axis=0).mask
      #print ma_average(Yz, axis=0), "mask=", ma_average(Yz, axis=0).mask
      #print sp.stats.pearsonr(ma_average(Xz, axis=0), ma_average(Yz, axis=0))
      #quit()

      # need +epsilon to avoid all zeros
      lsaTable[ti] = [i, j, Smax, Sl, Su, Xs, Ys, Al, Xs-Ys, lsaP, PCC, P_PCC, SPCC, P_SPCC, D_SPCC, SCC, P_SCC, SSCC, P_SSCC, D_SSCC]
      ti += 1
      del LSA_result
      #print "finalizing..."

  #print lsaTable
  #print "qvalues ..."
  #pvalues = np.ma.masked_invalid(pvalues)
  #pccpvalues = np.ma.masked_invalid(pccpvalues)
  #print "pvalues", pvalues, np.isnan(np.sum(pvalues))
  #print "pccpvalues", pccpvalues, np.isnan(np.sum(pccpvalues))
  #print "lsaP"
  #qvalues = storeyQvalue( pvalues )
  print >>sys.stderr, "LS Qvalues..."
  qvalues = qvalue_func( pvalues )
  #print "pccP"
  #pccqvalues = storeyQvalue( pccpvalues )
  print >>sys.stderr, "PCC Qvalues..."
  pccqvalues = qvalue_func( pccpvalues )
  #print "sccP"
  #sccqvalues = storeyQvalue( sccpvalues )
  print >>sys.stderr, "SCC Qvalues..."
  sccqvalues = qvalue_func( sccpvalues )
  #print "spccP"
  #spccqvalues = storeyQvalue( spccpvalues )
  print >>sys.stderr, "SPCC Qvalues..."
  spccqvalues = qvalue_func( spccpvalues )
  #print "ssccP"
  #ssccqvalues = storeyQvalue( ssccpvalues )
  print >>sys.stderr, "SSCC Qvalues..."
  ssccqvalues = qvalue_func( ssccpvalues )

  for k in xrange(0, len(qvalues)):
    lsaTable[k].append( qvalues[k] )
    lsaTable[k].append( pccqvalues[k] )
    lsaTable[k].append( spccqvalues[k] )
    lsaTable[k].append( sccqvalues[k] )
    lsaTable[k].append( ssccqvalues[k] )
    
  #print lsaTable
  for row in lsaTable:
    print >>resultFile, "\t".join(['%s']*len(col_labels)) % \
      tuple( [firstFactorLabels[row[0]], secondFactorLabels[row[1]]]  
          + ["%f" % np.round(v, decimals=disp_decimal) if isinstance(v, float) else v for v in row[2:]]
          + [row[0]+1, row[1]+1] )

### Liquid Association Section ###
def applyLA(inputData, scoutVars, factorLabels, bootCI=.95, bootNum=1000, minOccur=.50, pvalueMethod=1000,\
    fTransform=simpleAverage, zNormalize=noZeroNormalize, resultFile=None):

  col_labels = ['X','Y','Z','LA','lowCI','upCI','P','Q','Xi','Yi','Zi']
  print >>resultFile,  "\t".join(col_labels)

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
  for i in xrange(0, scoutNum):
    Xi = scoutVars[i][0] - 1
    Yi = scoutVars[i][1] - 1
    #print Xi, Yi
    Xo = np.ma.masked_invalid(inputData[Xi], copy=True) #need to convert to masked array with na's, not F-normalized
    Yo = np.ma.masked_invalid(inputData[Yi], copy=True) #need to convert to masked array with na's, not F-normalized
    for j in xrange(0, inputFactorNum):
      Zi = j
      if Xi == Yi or Xi == Zi or Zi == Yi:
        continue   #ignore invalid entries
      if cp[Xi,Yi,Zi] or cp[Xi,Zi,Yi] or cp[Zi,Xi,Yi] or cp[Zi,Yi,Xi] or cp[Yi,Xi,Zi] or cp[Yi,Zi,Xi]:
        continue   #ignore redundant entries
      cp[Xi,Yi,Zi]=True
      Zo = np.ma.masked_invalid(inputData[Zi], copy=True)    # need to convert to masked array with na's, not F-normalized
      Xo_minOccur = np.sum(np.logical_or(np.isnan(ma_average(Xo)), ma_average(Xo)==0))/float(timespots) < minOccur
      Yo_minOccur = np.sum(np.logical_or(np.isnan(ma_average(Yo)), ma_average(Yo)==0))/float(timespots) < minOccur
      Zo_minOccur = np.sum(np.logical_or(np.isnan(ma_average(Yo)), ma_average(Yo)==0))/float(timespots) < minOccur
      if Xo_minOccur or Yo_minOccur or Zo_minOccur:
        continue						 # minOccur not satisfied for one of the factors
      if np.all(Xo.mask) or np.all(Yo.mask) or np.all(Zo.mask):  # not any unmasked value in Xz or Yz, all nan in input, continue
        continue
      LA_score = singleLA(Xo, Yo, Zo, fTransform, zNormalize)                          # do LA computation

      #This Part Must Follow Static Calculation Part
      #Otherwise Yz may be changed, now it is copied
      #np.ma.array(copy=True) to copy, otherwise is only reference
      if pvalueMethod >= 0 or pvalueMethod < 0:
        Xp = np.ma.array(Xo,copy=True)
        Yp = np.ma.array(Yo,copy=True)
        Zp = np.ma.array(Zo,copy=True)
        laP = LApermuPvalue(Xp, Yp, Zp, pvalueMethod, np.abs(LA_score), fTransform, zNormalize)          # do Permutation Test
      else: #reserved 
        print "this should not be receahed"
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
  #qvalues = storeyQvalue( pvalues )
  qvalues = qvalue_func( pvalues )
  for k in xrange(0, len(qvalues)):
    laTable[k] = laTable[k] + [ qvalues[k], laTable[k][0]+1, laTable[k][1]+1, laTable[k][2]+1 ]

  #print laTable
  for row in laTable:
    print >>resultFile, "\t".join( ['%s']*len(col_labels) ) % \
      tuple( [factorLabels[row[0]], factorLabels[row[1]], factorLabels[row[2]]] \
          + ["%f" % np.round(v, decimals=disp_decimal) if isinstance(v, float) else v for v in row[3:]] )

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
    Xb = np.ma.array([ sample_wr(series1[:,j], series1.shape[0]) for j in xrange(0,series1.shape[1]) ]).T
    Yb = np.ma.array([ sample_wr(series2[:,j], series2.shape[0]) for j in xrange(0,series2.shape[1]) ]).T
    Zb = np.ma.array([ sample_wr(series3[:,j], series3.shape[0]) for j in xrange(0,series3.shape[1]) ]).T
    BS_set[i] = compcore.calc_LA(Xb, Yb, Zb)
  BS_set.sort()                                 #from smallest to largest
  BS_mean = np.mean(BS_set)
  return ( BS_mean, BS_set[np.floor(bootNum*a1)-1], BS_set[np.ceil(bootNum*a2)-1] )

def LApermuPvalue(series1, series2, series3, pvalueMethod, LA_score, fTransform, zNormalize):
  PP_set = np.zeros(pvalueMethod, dtype='float')
  X = zNormalize(fTransform(series1))
  Y = zNormalize(fTransform(series2))
  Z = np.ma.array(series3)                                               #use = only assigns reference, must use a constructor
  for i in xrange(0, pvalueMethod):
    np.random.shuffle(Z.T)
    PP_set[i] = compcore.calc_LA(X, Y, zNormalize(fTransform(Z)))
  if LA_score >= 0:
    P_two_tail = np.sum(np.abs(PP_set) >= LA_score)/np.float(pvalueMethod)
  else:
    P_two_tail = np.sum(-np.abs(PP_set) <= LA_score)/np.float(pvalueMethod)
  return P_two_tail

def test():
  """ self test script
  """
  np.seterr(all='raise')
  print >>sys.stderr, "###testing###"
  test_data = np.array( [ [ [ 3, 2, 0, 4], [1, 6, 0, 3], [6, 4, 0, 8], [np.nan, np.nan, np.nan, np.nan] ], 
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
  print >>sys.stderr, "---noZeroNormalize---"
  X=np.array( [1, 2, 3, np.nan, 0, 0], dtype='float')
  Xz=np.ma.masked_invalid(X)
  print >>sys.stderr, Xz, "->", noZeroNormalize(Xz)
  print >>sys.stderr, "---storeyQvalue---"
  pvalues = np.array([0.01, 0.2, 0.03, 0.4, 0.05, np.nan, 0.03, 0.4, 0.03, 0.3], dtype='float')
  print >>sys.stderr, "pvalues:", pvalues 
  print >>sys.stderr, "qvalues:", storeyQvalue(pvalues)
  pvalues = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.02, 0.03, 0.04, 0.03, 0.03], dtype='float')
  print >>sys.stderr, "pvalues:", pvalues 
  print >>sys.stderr, "qvalues:", storeyQvalue(pvalues)
  print >>sys.stderr, "---singleLSA---"
  print >>sys.stderr, "masked_data[0]:", masked_data[0]
  sa0 = simpleAverage(masked_data[0])
  print >>sys.stderr, "simpleAverage of masked_data[0]:", sa0, sa0.mask
  print >>sys.stderr, "input data:", noZeroNormalize(simpleAverage(masked_data[0]))
  print >>sys.stderr, "masked_data[1]:", masked_data[1]
  sa1 = simpleAverage(masked_data[1])
  print >>sys.stderr, "simpleAverage of masked_data[1]:", sa1, sa1.mask
  print >>sys.stderr, "input data:", noZeroNormalize(simpleAverage(masked_data[1]))
  lsar=singleLSA(masked_data[0], masked_data[1], delayLimit=1, fTransform=simpleAverage, zNormalize=noZeroNormalize, keepTrace=True)
  print >>sys.stderr, "lsar.score=", lsar.score 
  Al=len(lsar.trace)
  print >>sys.stderr, "lsar.align=",(lsar.trace[Al-1][0], lsar.trace[Al-1][1], len(lsar.trace)) 
  print >>sys.stderr, "---bootstrapCI---"
  print >>sys.stderr, "Bset=", bootstrapCI(masked_data[0], masked_data[1], lsar.score, 1, .95, 2, simpleAverage, noZeroNormalize)
  print >>sys.stderr, "---permuPvalue---"
  print >>sys.stderr, "P=", permuPvalue(masked_data[1], masked_data[0], 1, 2, np.abs(lsar.score), simpleAverage, noZeroNormalize)
  print >>sys.stderr, "---PCC---"
  (nPCC, nP_PCC) = sp.stats.pearsonr(np.mean(np.nan_to_num(test_data[0]), axis=0), np.mean(np.nan_to_num(test_data[1]), axis=0))
  oPCC = sp.corrcoef( np.mean(np.nan_to_num(test_data[0]),axis=0), np.mean(np.nan_to_num(test_data[1]),axis=0) )[0,1]
  otcdf = sp.stats.distributions.t.cdf(oPCC*np.sqrt((test_tN-2)/np.float(1.000000001-oPCC**2)), (test_tN-2))
  oP_PCC = .5 + np.sign(oPCC)*(.5 - otcdf) #addhoc for perfect correlation
  print >>sys.stderr, "nPCC", "nP_PCC", "oPCC", "oP_PCC", "otcdf"
  print >>sys.stderr, nPCC, nP_PCC, oPCC, oP_PCC, otcdf
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, clean_data, True, 1, .95, 10, 10, sdAverage, noZeroNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, clean_data, True, 1, .95, 10, 10, simpleAverage, noZeroNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, clean_data, True, 1, .95, 10, 10, simpleMedian, noZeroNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, clean_data, True, 1, .95, 10, 10, madMedian, noZeroNormalize)
  #print >>sys.stderr, "---theoPvalue--- d=3, xi=(0 to 100)*100, x=(xi/100)*x_var"
  #D3_table = theoPvalue(Rmax=10, Dmax=3, precision=.00001, x_decimal=my_decimal)   #let's produce 2 tail-ed p-value
  #D2_table = theoPvalue(Rmax=10, Dmax=2, precision=.00001, x_decimal=my_decimal)   #let's produce 2 tail-ed p-value
  #D1_table = theoPvalue(Rmax=10, Dmax=1, precision=.00001, x_decimal=my_decimal)   #let's produce 2 tail-ed p-value
  #D0_table = theoPvalue(Rmax=10, Dmax=0, precision=.00001, x_decimal=my_decimal)   #let's produce 2 tail-ed p-value
  #numpy.arange(2,5.1,0.2) x=R/sqrt(N)
  #for r in np.arange(0, R, 0.01):
  #print readPvalue(D1_table, R, N, x_sd=1, M=1, alpha=1, beta=1, x_decimal=3):
  #print >>open("D0.txt",'w'), "\n".join([ "%s\t%s" % (str(v/float(10**my_decimal)), str(T_table[v])) for v in T_table.keys() ])
  #print >>open("D1.txt",'w'), "\n".join([ "%s\t%s" % (str(v/float(10**my_decimal)), str(T_table[v])) for v in T_table.keys() ])
  #print >>open("D2.txt",'w'), "\n".join([ "%s\t%s" % (str(v/float(10**my_decimal)), str(T_table[v])) for v in T_table.keys() ])
  #print >>open("D3.txt",'w'), "\n".join([ "%s\t%s" % (str(v/float(10**my_decimal)), str(T_table[v])) for v in T_table.keys() ])
  #print >>sys.stderr, P=readPvalue(T_table, .1876, 140, x_var=1, x_decimal=my_decimal) # read two-tailed  
  #print >>sys.stderr, "theoP=", P

if __name__=="__main__":
  print "hello world!"
  test()
  exit(0)


#############  FUNCTIONS NO IN USE AFTER THIS LINE ####################

