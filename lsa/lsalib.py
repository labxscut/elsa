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
			zNormalize(func):	nomalizing function

		Return:
			one single LSA result
      new implemented similarity alignment using external C++ routine in compcore
    
  """

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

def bootstrapCI(series1, series2, delayLimit, bootCI, bootNum, fTransform, zNormalize):
  """	do bootstrap CI estimation

		Args:
			series1(np.array): 	sequence data of Seq X
			series2(np.array): 	sequence data of Seq Y
			delayLimit(int): 	maximum time unit of delayed response allowed	
      bootCI(float):  confidence interval size
      bootNum(int): number of bootstraps
			fTransform(func):	replicate summarizing function
			zNormalize(func):	nomalizing function

    Return:
      Confidence Interval

	"""
	
  ###print "------Bootstrapping------"
  lsad = compcore.LSA_Data()
  BS_set = np.zeros(bootNum, dtype='float')
  for i in range(0, bootNum):
    Xb = np.array([ sample_wr(series1[:,j], series1.shape[0]) for j in xrange(0,series1.shape[1]) ]).T
    Yb = np.array([ sample_wr(series2[:,j], series2.shape[0]) for j in xrange(0,series2.shape[1]) ]).T
    lsad.assign( delayLimit, zNormalize(fTransform(Xb)), zNormalize(fTransform(Yb)) )
    BS_set[i] = compcore.DP_lsa(lsad, False).score
  BS_set.sort()
  #print "BS_set=", BS_set
  return ( BS_set[np.floor(bootNum*(1-bootCI)/2.0)], BS_set[np.floor(bootNum*bootCI/2.0)] )
	
def permuPvalue(series1, series2, delayLimit, permuNum, Smax, fTransform, zNormalize):
  """ do permutation Test

		Args:
			series1(np.array): 	sequence data of Seq X
			series2(np.array): 	sequence data of Seq Y
			delayLimit(int): 	maximum time unit of delayed response allowed	
      permuNum(int): number of permutations
      Smax(int): maximum LSA
			fTransform(func):	replicate summarizing function
			zNormalize(func):	nomalizing function

    Return:
      p-value

	"""
	
  ###print "-------permutation------"
  lsad = compcore.LSA_Data()
  PP_set = np.zeros(permuNum, dtype='float')
  Xz = zNormalize(fTransform(series1))
  Y = np.array(series2)                                               #use = only assigns reference, must use a constructor
  #print "series2=", series2, fTransform(series2), zNormalize(fTransform(series2))
  for i in xrange(0, permuNum):
    np.random.shuffle(Y.T)
    #print "Y=", Y, fTransform(Y), zNormalize(fTransform(Y))
    lsad.assign( delayLimit, Xz, zNormalize(fTransform(Y)) )
    PP_set[i] = np.abs(compcore.DP_lsa(lsad, False).score)
  #print "series2=", series2, fTransform(series2), zNormalize(fTransform(series2))
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
    mse = np.zeros(100, dtype='float')
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
  """

  return np.average(tseries, axis=0)


def	sdAverage(tseries):
  """	sd weighted averaging 

    Args:
      tseries(np.array):  one time series with replicates, each row is a replicate

    Reterns:
      one row with replicates sd weighted averaged
  """

  return np.average(tseries, axis=0)/np.std(tseries,axis=0) 


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
  return np.array( [ v_cum[values[i]] for i in xrange(0, len(values)) ], dtype='float' )


def	scoreNormalize(tseries):
  """	score normalizing

    Args:
      tseries(np.array):  one time series with no replicates

    Reterns:
      score normalized tseries
  """
  p_ranks = tied_rank(tseries)
  return sp.stats.distributions.norm.ppf( p_ranks/(len(p_ranks)+1) )

def fillMissing(tseries, method):
  """ fill missing data

    Args:
      tseries(np.array):  one time series with no replicates
      method(str):  filling method, choices ['none', 'zero', 'linear', 'slinear', 'nearest', 'quadratic', 'cubic'] 

    Reterns:
      tseries with missing data filled
  """
  
  if method == 'none':
    return np.nan_to_num(tseries)
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
    
def applyAnalysis( cleanData, delayLimit=3, bootCI=.95, bootNum=1000, permuNum=1000, fTransform=simpleAverage, zNormalize=scoreNormalize ):
  """ calculate pairwise LS scores and p-values

    	Args:
    		cleanData(np.array): 	numpy data array with correct format factor_num x timespot_num x replicte_num, no missing values
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

  #print factorNum, repNum, spotNum, lsaTable, pvalues

  ti = 0
  for i in xrange(0, factorNum-1):
    for j in xrange(i+1, factorNum):
      LSA_result = singleLSA(cleanData[i], cleanData[j], delayLimit, fTransform, zNormalize, True)                          # do LSA computation
      Smax = LSA_result.score                                                                                         # get Smax
      (Sl, Su) = bootstrapCI(cleanData[i], cleanData[j], delayLimit, bootCI, bootNum, fTransform, zNormalize)         # do Bootstrap CI
      permuP = permuPvalue(cleanData[i], cleanData[j], delayLimit, permuNum, np.abs(Smax), fTransform, zNormalize)    # do Permutation Test
      pvalues[ti] = permuP
      Al = len(LSA_result.trace)
      (Xs, Ys, Al) = (LSA_result.trace[Al-1][0], LSA_result.trace[Al-1][1], len(LSA_result.trace))
      #print cleanData[i], cleanData[j], np.mean(cleanData[i], axis=0)+np.finfo(np.double).eps, np.mean(cleanData[j], axis=0)+np.finfo(np.double).eps
      (PCC, P_PCC) = sp.stats.pearsonr(np.mean(cleanData[i], axis=0), np.mean(cleanData[j], axis=0))# +epsilon to avoid all zeros
      #PCC = sp.corrcoef( cleanData[i], cleanData[j] )[0,1]
      #tcdf = sp.stats.distributions.t.cdf( PCC*np.sqrt((spotNum-1)/np.float(1.000000001-PCC**2)), (spotNum-1))
      #P_PCC = .5 + np.sign(PCC)*(.5 - tcdf )                                                                             #addhoc for perfect correlation
      lsaTable[ti] = [i, j, Smax, Sl, Su, Xs, Ys, Al, Xs-Ys, permuP, PCC, P_PCC]
      ti += 1

  #print "pvalues", pvalues
  qvalues = storeyQvalue( pvalues )
  for i in xrange(0, len(qvalues)):
    lsaTable[i].append( qvalues[i] )

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
      clean_data[i][j] = fillMissing(test_data[i][j], 'zero')
  print >>sys.stderr, "clean_data", clean_data
  print >>sys.stderr, "---tied-rank---"
  print >>sys.stderr, tied_rank(clean_data[0][0])
  print >>sys.stderr, "---scoreNormalize---"
  print >>sys.stderr, scoreNormalize(clean_data[0][0])
  print >>sys.stderr, "---simpleAverage---" 
  print >>sys.stderr, simpleAverage(clean_data[0])
  print >>sys.stderr, "---sdAverage---"
  print >>sys.stderr, sdAverage(clean_data[0])
  print >>sys.stderr, "---storeyQvalue---"
  pvalues = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.02, 0.03, 0.04, 0.03, 0.03], dtype='float')
  print >>sys.stderr, "pvalues:", pvalues 
  print >>sys.stderr, storeyQvalue(pvalues)
  print >>sys.stderr, "---singleLSA---"
  lsar = singleLSA(clean_data[0], clean_data[1], delayLimit=1, fTransform=simpleAverage, zNormalize=scoreNormalize, keepTrace=True)
  print >>sys.stderr, lsar.score
  print >>sys.stderr, "---bootstrapCI---"
  print >>sys.stderr, bootstrapCI(clean_data[0], clean_data[1], 1, .99, 10, simpleAverage, scoreNormalize)
  print >>sys.stderr, "---permuPvalue---"
  print >>sys.stderr, "p-value:", permuPvalue(clean_data[1], clean_data[0], 1, 10, np.abs(lsar.score), simpleAverage, scoreNormalize)
  print >>sys.stderr, "---applyAnalysis---"
  print >>sys.stderr, applyAnalysis(clean_data, 1, .99, 10, 10, simpleAverage, scoreNormalize)
  

if __name__=="__main__":
  print "hello world!"
  test()
  exit(0)


#############  FUNCTIONS NO IN USE AFTER THIS LINE ####################

##########################################            
#### NO NEED to look beyond this line ####
##########################################

######################################
#linearFillRow						##	
#TODO: 								##
######################################
#def linearFillRow( csvReader, csvWriter, repnum=1, placeHolder="na", mode="zero", skiprows=1):
#    """linearly fill the missing value by in/extrapolation, serie data in column.
#
#    csvReader - csv file descriptor of the raw data file
#    csvWriter - csv file descriptor of a file for temporary storage of data
#    repnum - how many replicates, default 1
#    placeHolder - place holder in data file for missing value, default: na
#    mode - in/extrapolation to be used, default: intropolation only
#    skiprows - rows to skip before data
#    """
#    idx=0
#    for row in csvReader:
#        if idx < skiprows:
#        	#skip the labels
#            csvWriter.writerow(row)
#            idx = idx + 1
#        else:
#        	#retrieve the data
#        	colnum = len(row) - 1
#        	replicates = np.zeros((repnum, colnum/repnum), dtype='float').tolist()
#        	for col in range(1, colnum):
#        		#get independent replicate
#        		replicates[col/repnum][col-col/repnum] = row[col]
#        	for rep in replicates:
#        		#fix them one by one
#            	for i in range(0, repnum):
#                	j = i
#                	if rep[i] == "na" or rep[i] == '':	
#                    	while j < len(rep):
#                        	if rep[j] != "na" and rep[j] != '':
#                            	break	#find na's
#                        	j=j+1
#                	if i == 0: # na from start
#                    	if mode == "in/ex" and j <= len(rep)-2 and (rep[j+1] != "na" and rep[j+1] != ''): # inex mode and extropolate possible
#                        	for k in range(i,j):
#                            	rep[k]="%4.f" % (float(rep[j])+(float(rep[j])-float(rep[j+1]))*(j-k))
#                    	else:
#                        	for k in range(i,j):
#                            	rep[k]="0"
#                	elif j == len(rep): # na at end
#                    	if mode == "in/ex" and i >= 2: # inex mode and extropolate possible
#                        	for k in range(i,j):
#                            	rep[k]="%.4f" % (float(rep[i-1])+(float(rep[i-1])-float(rep[i-2]))*(k-(i-1)))
#                    	else:
#                        	for k in range(i,j):
#                            	rep[k]="0"
#                	else: # intropolate or zero
#                    	for k in range(i,j):
#                        	if mode != "zero":
#                            	rep[k]="%.4f" % (float(rep[i-1])+(float(rep[j])-float(rep[i-1]))/(j-i+1)*(k-i+1))
#                        	else:
#                            	rep[k]="0"
#                	i = j
#            for col in range(1, colnum):
#            	#reconstruct the data row
#            	row[col] = replicates[col/repnum][col-col/repnum]
#            csvWriter.writerow(row) 


###############################
# sigTestMaxDelay
# return a list-of-list, where i-th list = [ s1, s2, L-score, X, Y, length, P/N, P-value, corr, corPval
############################### 

#def newSigTestMaxDelay( dataX, maxDelay=3, permu=1000 ):   
#
#    
#    lsMatrix = []
#    counter = 0
#    for i in range( dataX.shape[1]-1 ):
#        for j in range( i+1, dataX.shape[1] ):
#            pairLs = newLocalSimilarity( dataX[:,i], dataX[:,j], maxDelay, False )
#            pValue = newSigTest( dataX[:,i], dataX[:,j], pairLs[0], maxDelay, permu, False )
#            corr = corrcoef( dataX[:,i], dataX[:,j] )[0,1] 
#            tcdf = t.cdf( corr*sqrt((dataX.shape[0]-1)/(float64)(1.000000001-corr**2)), (dataX.shape[0]-1)) 
#            #adhot approach for dealing with perfect correlation
#            corpVal = .5 + sign(corr)*(.5 - tcdf )  
#            pairLs[0] = "%.4f" % ((pairLs[0]/(float)(dataX.shape[0]))) # normalize and add sign
#            pValue = "%.5f" % pValue
#            corr = "%.4f" % corr
#            corpVal = "%.5f" % corpVal
#            lsMatrix.append( [i+1, j+1] + pairLs[0:4] + [pairLs[1]-pairLs[2]] + [pValue, corr, corpVal] )
#            counter = counter + 1
#    return lsMatrix


#################################
# normalTransform performs rank normal score transform (Ker-chau Li, PNAS 2002)
# rawMatrix is an original data matrix, with each column as factor and each row as a time point
# rankScoreMatrix(i,:)=inverse normal cumulative distribution of rank(rawMatrix(i,:))/(rowNum+1)
#################################

#def normalTransform( rawMatrix ):
#    """normal transformation of the raw data
#
#        rawMatrix - matrix of the raw data
#    """
#    rankScoreMatrix = -ones(rawMatrix.shape, dtype=float64)
#    for i in range(rawMatrix.shape[1]):
#        rankScoreMatrix[:,i] = norm.ppf( stats.rankdata(rawMatrix[:,i])/(rawMatrix.shape[0]+1) )
#    return rankScoreMatrix
#
#def rowNormalTransform( rawMatrix ):
#    rankScoreMatrix = -ones(rawMatrix.shape, dtype=float64)
#    for i in range(rawMatrix.shape[0]):
#      rankScoreMatrix[i,:] = norm.ppf( stats.rankdata(rawMatrix[i,:])/(rawMatrix.shape[1]+1) )
#    return rankScoreMatrix

###############################
# Z-normalization
############################### 

#def dataNormalize( dataX ): # assume missing value has already been dealt
#    """normal(Z) transformation of raw data
#        
#        dataX - matrix of sequence data, raw data
#    """
#    if rank(dataX) == 1:
#        dataX = (dataX - mean(dataX))/std(dataX)
#    if rank(dataX) == 2:
#        for i in range(dataX.shape[1]):
#            dataX[:,i] = (dataX[:,i] - (dataX[:,i]).mean())/(dataX[:,i]).std()
#    return dataX
#
#def rowDataNormalize( dataX ):
#    if rank(dataX) == 1:
#        dataX = (dataX - mean(dataX))/std(dataX)
#    if rank(dataX) == 2:
#        for i in range(dataX.shape[0]):
#            dataX[i,:] = (dataX[i,:] - (dataX[i,:]).mean())/(dataX[i,:]).std()
#    return dataX
    
###############################
# localSimilarity
# return a list = [ L-score, start-X, start-Y, length, Ahead/Lag ]
############################### 

#def newLocalSimilarity( lsTS1, lsTS2, maxDelay=3, scale=False ):
#    """do local simularity alignment and return a list in following format:
#        [ Local Similarity(LS) Score, Seq X's Align Start Position, 
#            Seq Y's Align Start Position, Length of Aligned Potion, 
#            Positive/Negtive Correlated ]
#        new implemented similarity alignment using external C++ routine
#
#    lsTS1 - sequence data of Seq X
#    lsTS2 - sequence data of Seq Y
#    maxDelay - maximum time unit of delayed response allowed
#    scale - scale the data before alignment, default: False
#    """
#    lsad=compcore.LSA_Data(maxDelay,lsTS1,lsTS2)
#    lsar=compcore.DP_lsa(lsad)
#    len=lsar.trace.size()
#    if lsar.score > 0:
#        posOrNeg = 1
#    elif lsar.score < 0:
#        posOrNeg = -1
#    else:
#        posOrNeg = 0
#    #for i in range(0, len):
#    #  print lsar.trace[i][0], "\t", lsar.trace[i][1]
#    #quit()
#    return [lsar.score, lsar.trace[len-1][0], lsar.trace[len-1][1], len, posOrNeg]

###############################
# sigTest
# return a float = significance(P-value) based on a permutation test
############################### 


#def newSigTest( lsTS1, lsTS2, refScore=0, maxDelay=3, permu=1000, scale=False ):
#    """significance test by permutation, return a P-value for the LS score.
#    lsTS1 - sequence data of Seq X
#    lsTS2 - sequence data of Seq Y
#    refScore - score from local similarity alignment
#    maxDelay - maximum time unit of delayed response allowed
#    permu - number of permutation time
#    scale - scale the data before alignment, default: False
#
#    using new C++ routine 
#    """
#    lsad=compcore.LSA_Data(maxDelay,lsTS1,lsTS2)
#    lsar=compcore.DP_lsa(lsad)
#    lsa_dpr=compcore.LSA_PT(lsad, lsar, permu, compcore.LSA_test)
#    return lsa_dpr.pvalue
    
 
#############################################
### BELOW, Functions not in use currently ###
### MOSTLY, old python based functions    ###
#############################################

#def linearFillColumn( csvReader, csvWriter, placeHolder="na", mode="zero", skiprows=1): 
#    """linearly fill the missing value by in/extrapolation, serie data in column.
#
#    csvReader - csv file descriptor of the raw data file
#    csvWriter - csv file descriptor of a file for temporary storage of data
#    placeHolder - place holder in data file for missing value, default: na
#    mode - in/extrapolation to be used, default: intropolation only
#    skiprows - rows to skip before data
#    """
#    table=[]
#    idx=0
#    for row in csvReader:
#        if idx < skiprows:
#            csvWriter.writerow(row)
#            idx = idx + 1
#            continue
#        else:
#            table.append(row)
#    # exception handler: table empty?
#    for l in range(1, len(table[0])):
#        for i in range(0,len(table)):
#            j = i
#            if table[i][l] == "na" or table[i][1] == '':
#                while j < len(table): 
#                    if table[j][l] != "na" and table[j][1] != '': 
#                        break
#                    j = j + 1
#            if i == 0: # na from start
#                if mode == "in/ex" and j <= len(table)-2 and (table[j+1][l] != "na" and table[j+1][1] != ''): # inex mode and extropolate possible
#                    for k in range(i,j):
#                        table[k][l]="%.4f" % (float(table[j][l])+(float(table[j][l])-float(table[j+1][l]))*(j-k))
#                else:
#                    for k in range(i,j):
#                        table[k][l]="0"
#            elif j == len(table): # na from end
#                if mode == "in/ex" and i >= 2: # inex mode and extropolate possible
#                    for k in range(i,j):
#                        table[k][l]="%.4f" % (float(table[i-1][l])+(float(table[i-1][l])-float(table[i-2][l]))*(k-(i-1)))
#                else:
#                    for k in range(i,j):
#                        table[k][l]="0"
#            else: # intropolate
#                for k in range(i,j):
#                    if mode != "zero":
#                        table[k][l]="%.4f" % (float(table[i-1][l])+(float(table[j][l])-float(table[i-1][l]))/(j-i+1)*(k-i+1))
#                    else:
#                        table[k][l]="0"
#            i = j
#    csvWriter.writerows(table)
#    
#def sigTestMaxDelay( dataX, maxDelay=3, permu=1000 ):   
#    """calculate pair to pair LS and P-value and return a lsa table.
#    each row in the table is a list in following format:
#    [ Seq X's Idx, Seq Y's Idx, LS Score, X's Start Position, 
#        Y's Start Position, Alignment Length, Delay in Time Unit,
#        P-value, Pearson' Correlation, Correlation P-value ]
#
#    dataX - transformed data matrix
#    maxDelay - maximum time unit of delayed response allowed
#    permu - number of permutation time
#    """
#    
#    lsMatrix = []
#    counter = 0
#    for i in range( dataX.shape[1]-1 ):
#        for j in range( i+1, dataX.shape[1] ):
#            pairLs = localSimilarity( dataX[:,i], dataX[:,j], maxDelay, False )
#            pValue = sigTest( dataX[:,i], dataX[:,j], pairLs[0], maxDelay, permu, False )
#            corr = corrcoef( dataX[:,i], dataX[:,j] )[0,1] 
#            # or use sd.correlation(X,Y), see discussion online by search "correlation coefficient numpy"
#            #corpVal = .5 + sign(corr)*(.5 - t.cdf( corr*sqrt((dataX.shape[0]-1)/(float)(1-corr**2)), dataX .shape[0]-1))  
#            tcdf = t.cdf( corr*sqrt((dataX.shape[0]-1)/(float64)(1.000000001-corr**2)), (dataX.shape[0]-1))                 
#            #adhot approach for dealing with perfect correlation
#            corpVal = .5 + sign(corr)*(.5 - tcdf )
#            #corpVal = .5 + sign(corr)*(.5 - t.cdf( corr*sqrt((dataX.shape[0]-1)/(float)(1-corr**2)), \
#            #  (float)(dataX.shape[0]-1)))  
#            pairLs[0] = "%.4f" % (pairLs[0]/(float)(dataX.shape[0])*int(pairLs[4])) # normalize and already signed
#            pValue = "%.5f" % pValue
#            corr = "%.4f" % corr
#            corpVal = "%.5f" % corpVal
#            lsMatrix.append( [i+1, j+1] + pairLs[0:4] + [pairLs[1]-pairLs[2]] + [pValue, corr, corpVal] )
#            counter = counter + 1
#    return lsMatrix
#
#def localSimilarity( lsTS1, lsTS2, maxDelay=3, scale=False ):
#    """do local simularity alignment and return a list in following format:
#        [ Local Similarity(LS) Score, Seq X's Align Start Position, 
#            Seq Y's Align Start Position, Length of Aligned Potion, 
#            Positive/Negtive Correlated ]
#
#    lsTS1 - sequence data of Seq X
#    lsTS2 - sequence data of Seq Y
#    maxDelay - maximum time unit of delayed response allowed
#    scale - scale the data before alignment, default: False
#    """
#    #if scale == True:
#    #   lsTS1 = (lsTS1-mean(lsTS1))/std(lsTS1)
#    #   lsTS2 = (lsTS2-mean(lsTS2))/std(lsTS2)
#    posScoreMatrix = zeros( (lsTS1.shape[0]+1,lsTS2.shape[0]+1), dtype=float64 )
#    negScoreMatrix = zeros( (lsTS1.shape[0]+1,lsTS2.shape[0]+1), dtype=float64 )
#    scoreMax = 0.
#    posOrNeg = 0 # 1 for positive, 0 for negative, -1 for unknown 
#    startX = 0
#    startY = 0
#    for i in range(1, lsTS1.shape[0]+1):
#        #for j in range(1, lsTS2.shape[0]+1):
#        for j in range(maximum(1, i-maxDelay), minimum(lsTS2.shape[0]+1, i+maxDelay+1)):
#            #if abs(i-j) > maxDelay:
#            #   continue
#            #else:
#            posScoreMatrix[i,j]=max(0, posScoreMatrix[i-1,j-1]+lsTS1[i-1]*lsTS2[j-1])
#            negScoreMatrix[i,j]=max(0, negScoreMatrix[i-1,j-1]-lsTS1[i-1]*lsTS2[j-1])
#            if posScoreMatrix[i,j] > scoreMax:
#                scoreMax = posScoreMatrix[i,j]
#                startX = i
#                startY = j
#                posOrNeg = 1
#            if negScoreMatrix[i,j] > scoreMax:
#                scoreMax = negScoreMatrix[i,j]
#                startX = i
#                startY = j
#                posOrNeg = -1
#    #thresh = .00001
#    length = 0
#    if posOrNeg == 1:
#        while posScoreMatrix[startX - length, startY - length] != 0:
#            length = length + 1
#    elif posOrNeg == -1:
#        while negScoreMatrix[startX - length, startY - length] != 0:
#            length = length + 1
#    #pdb.set_trace()
#    return [scoreMax, startX - length + 1, startY - length + 1, length, posOrNeg]
#    
#
#def sigTest( lsTS1, lsTS2, refScore=0, maxDelay=3, permu=1000, scale=False ):
#    """significance test by permutation, return a P-value for the LS score.
#    lsTS1 - sequence data of Seq X
#    lsTS2 - sequence data of Seq Y
#    refScore - score from local similarity alignment
#    maxDelay - maximum time unit of delayed response allowed
#    permu - number of permutation time
#    scale - scale the data before alignment, default: False
#    """
#    rTS1=zeros_like(lsTS1)+lsTS1
#    rTS2=zeros_like(lsTS2)+lsTS2
#    highScoreNum = 0
#    for i in range(permu):
#        shuffle(rTS1) 
#        shuffle(rTS2) 
#        scoreRandom = localSimilarity( rTS1, rTS2, maxDelay, scale )[0]
#        if scoreRandom >= refScore:
#            highScoreNum = highScoreNum + 1
#    pValue = highScoreNum / (float64)(permu)
#    return pValue

###############################
# temporalChangeSeq 
#
###############################

#def temporalChangeSeq( dataX ):
#    """calculate the temporal change between pairs of sequences
#
#        dataX - matrix of sequence data, raw data
#    """
#    data1 = zeros((dataX.shape[0], dataX.shape[1]-1), dtype=float64)
#    data2 = zeros(data1.shape, dtype=float64)
#    data1 = dataX[:,1:dataX.shape[1]] - dataX[:,0:dataX.shape[1]-1]
#    data2 = dataX[:,1:dataX.shape[1]] + dataX[:,0:dataX.shape[1]-1]
#    data2 = where(data2==0, 108, data2) # not ought to be 108, any value except 0 will work
#    return data1/data2
            
###############################
# 
#
############################### 

#def colSum20ne( dataX ):
#    """normalize data by column sum
#        
#        dataX - matrix of sequence data, either raw data or transformed
#    """
#    dataP = zeros(dataX.shape)
#    csum1 = dataX.sum(0) 
#    dataP = dataX / csum1
#    return dataP    
