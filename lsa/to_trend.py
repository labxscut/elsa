#!/usr/bin/env python

import argparse, time, sys
import numpy as np

try:
  #debug import
  import lsalib
  #np.seterr(all='raise')
except ImportError:
  #installed import
  from lsa import lsalib

def main():

  parser = argparse.ArgumentParser(\
      description="Tool Convert Profile Series to Trend Series")
  parser.add_argument("dataFile", metavar= "dataFile", \
      type=argparse.FileType('rU'), \
      help="the input data filename, use null for simulation")
  parser.add_argument("resultFile", metavar= "resultFile", \
      type=argparse.FileType('w'), \
      help="the output result file")
  parser.add_argument("-T", "--trendThresh", dest="trendThresh", default=0, \
      help="if specified must be a number, will generate trend series with the specified threshold")
  parser.add_argument("-S", "--spotNum", dest="spotNum", \
      help="time spot number", default=0)
  parser.add_argument("-N", "--normMethod", dest="normMethod", \
      help="normalization method", \
      choices=["none", "percentileZ", "percentile", "noZero"], default="none")
  arg_namespace = parser.parse_args()
  # parse arguments: 
  dataFile = vars(arg_namespace)['dataFile']
  resultFile = vars(arg_namespace)['resultFile']
  trend_threshold = float(vars(arg_namespace)['trendThresh'])
  spotNum = int(vars(arg_namespace)['spotNum'])
  normMethod = vars(arg_namespace)['normMethod']

  realData=np.genfromtxt( dataFile, comments='#', delimiter='\t', \
      missing_values=['na',''], filling_values=np.nan, \
      usecols=range(1,spotNum+1) )
  dataFile.seek(0)  #rewind
  factorLabels=list(np.genfromtxt( dataFile, comments='#', delimiter='\t', \
      usecols=xrange(0,1), dtype='string' ))

  start_time = time.time()

  #print factorLabels
  #print realData.shape
  rows, cols = realData.shape
  bootNum = 100
  #bootNum = 1
  #print rows, cols
  #print realData

  if trend_threshold==0:
    sigma_square = 1.25
  else:
    P = lsalib.calc_tmatrix(bootNum, trend_threshold)
    sigma_square = lsalib.calc_markov_var(P)
    # return parameters P=(a,b,c,d,t)
    # w, vl, vr = lsalib.calc_eigen(P)
    # rightEigenVec 3 by 1, leftEigenVec 1 by 3, r1 = 1, 
    # l1 = stationary distribution, lambda1=1
    # sigma_square = lsalib.calc_sigma_square(w, vl, vr)
  print >>sys.stderr, "sigma=", sigma_square

  trendData = []
  for j in range(0,rows):
    # trendify
    if normMethod == "none": # if already normalized
      OSeries = realData[j,:]
    elif normMethod == "percentileZ":
      OSeries = lsalib.percentileZNormalize(np.ma.masked_invalid(realData[j,:]))
    elif normMethod == "percentile":
      OSeries = lsalib.percentileNormalize(np.ma.masked_invalid(realData[j,:]))
    elif normMethod == "noZeroNormalize":
      OSeries = lsallib.noZeroNormalize(np.ma.masked_invalid(realData[j,:]))
     # print realData[j,:]
    TSeries = lsalib.ji_calc_trend(OSeries, cols, trend_threshold)
    trendData.append(TSeries)

  #print len(trendData), len(trendData[1])
  print >>resultFile, '\t'.join(['#V=%.2f' % sigma_square] \
      + [str(v) for v in range(1,cols)])
  for j in range(0, rows):
    print >>resultFile, '\t'.join([factorLabels[j]] + \
      ["na" if np.isnan(v) else "%d" % v for v in trendData[j]])

  end_time = time.time()
  elapse_time = end_time - start_time

  print >>sys.stderr, "finished in %d seconds" % elapse_time

if __name__=="__main__":
  main()
  exit(0)

#def ji_calc_trend(oSeries, lengthSeries, thresh):
#  #Liping Ji and Kian-Lee Tan, Bioinformatics 2005
#
#  #print oSeries
#  #print lengthSeries
#  tSeries = np.zeros(lengthSeries-1, dtype='float')
#
#  for i in xrange(0, lengthSeries-1):
#    if oSeries[i] == 0 and oSeries[i+1] > 0:
#      trend = 1
#    elif oSeries[i] == 0 and oSeries[i+1] < 0:
#      trend = -1
#    elif oSeries[i] == 0 and oSeries[i+1] == 0:
#      trend = 0
#    else:
#      trend = (oSeries[i+1]-oSeries[i])/np.abs(oSeries[i])
#
#    if np.isnan(trend):
#      tSeries[i] = np.nan
#    elif trend >= thresh:
#      tSeries[i] = 1
#    elif trend <= -thresh:
#      tSeries[i] = -1
#    else:
#      tSeries[i] = 0
#
#  #print tSeries
#  #exit()
#  return tSeries
#
#def calc_tmatrix(bootNum, trend_threshold, timeNum=1000, randomFunc=np.random.normal): 
#  # return a 3 by 3 transition matrix
#  Tm = np.zeros((bootNum,5)) #each row in order: a,b,c,d,t
#  for b in range(0, bootNum):
#    Tm[b,] = to_markov(trend_threshold, timeNum, rundomFunc)
#  print np.average(Tm, axis=1)
#  return np.average(Tm, axis=1)
#
#def calc_eigen(Tmatrix):
#  #rightEigenVec 3 by 1, leftEigenVec 1 by 3, r1 = 1, l1 = stationary distribution, lambda1=1
#  eigenValue, rightEigenVec, leftEigenVec = 0, 0, 0
#  pass
#
#def calc_sigma_square(eigenValue, rightEigenVec, leftEigenVec):
#  pass
#
#def to_markov(threshold, timeNum, randomFunc=np.random.normal): #t is threshold vector
#  N = timeNum
#  Xo = randomFunc(size=N+2)  #N+1 observations lead to N trends
#  Xt = []
#  for i in range(0,N+1):
#    if Xo[i] == 0 and Xo[i+1] > 0:
#      Xt.append(1)
#    elif Xo[i] == 0 and Xo[i+1] < 0:
#      Xt.append(-1)
#    elif Xo[i] == 0 and Xo[i+1] == 0:
#      Xt.append(0)
#    elif (Xo[i+1]-Xo[i])/np.abs(Xo[i]) >= t[j]: 
#      Xt.append(1)
#    elif (Xo[i+1]-Xo[i])/np.abs(Xo[i]) <= -t[j]:
#      Xt.append(-1)
#    else: 
#      Xt.append(0)
#  #print(Xo)
#  #print(Xt)
#  P = np.zeros((1,5)) 
#  for j in range(0,N):
#    #print(Xt[i,j])
#    if Xt[j]==1 and Xt[j+1]==1:
#      P[2] = P[2] +1
#    elif Xt[j]==1 and Xt[j+1]==-1:
#      P[3] = P[3] +1
#    elif Xt[j]==0 and Xt[j+1]==1:
#      P[4] = P[4] +1
#    else:
#      print "not possible"
#      pass
#    P[1] = np.sum(Xt[i] == 1) #di=1
#    P[1] = P[1]/(N+1)
#    P[2] = (P[2]/N)/(P[1])
#    P[3] = (P[3]/N)/(P[1])
#    P[4] = (P[4]/N)/(1-2*P[1])
#    P[5] = t[i]
#  return(P)

