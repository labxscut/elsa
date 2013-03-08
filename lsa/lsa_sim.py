#!/usr/bin/env python

import argparse, time, sys
import numpy as np

#kcut_min=100
#Rmax_min=10
my_decimal = 3    # preset x step size for P_table

try:
  #debug import
  import lsalib
  #np.seterr(all='raise')
except ImportError:
  #installed import
  from lsa import lsalib

def main():

  # define arguments: delayLimit, fillMethod, permuNum
  parser = argparse.ArgumentParser(description="LSA and LTA Simulation Tool")
  parser.add_argument("resultFile", metavar= "resultFile", \
      type=argparse.FileType('w'), help="the output result file")
  parser.add_argument("-S", "--simTimes", dest="simTimes", default=10000, type=int,\
      help="specify the times of simulation to run, default: 10000")
  parser.add_argument("-D", "--delayLimit", dest="delayLimit", default=0, type=int,\
      help="specify the maximum delay possible, default: 0")
  parser.add_argument("-L", "--lengthSeries", dest="lengthSeries", \
      default=50, type=int,\
      help="specify the length of series to generate, default: 50")
  parser.add_argument("-T", "--trendSeries", dest="trendSeries", default=None,\
      help="if specified must be a threshold number, will generate trend series, \
      with the specified threshold")
  parser.add_argument("-A", "--approxVar", dest="approxVar", default=1, type=float,\
      help="numeric>0, default=1, variance of partial sum variable")
  parser.add_argument("-a", "--alphaValue", dest="alphaValue", default=1, type=float,\
      help="1>=numeric>0, default=1, proportions of non-zeroes in X")
  parser.add_argument("-b", "--betaValue", dest="betaValue", default=1, type=float,\
      help="1>=numeric>0, default=1, proportions of non-zeroes in Y")
  parser.add_argument("-M", "--simMethod", dest="simMethod", default="idn,0,1",\
      help="idn,x,y: independent normal with mean x variance y")
  parser.add_argument("-n", "--nullDistribution", dest="nullDistribution", \
      default="yes", choices=['yes','no'],
      help="yes: randomize timpoints; no: do not randomize timepoints")
  parser.add_argument("-P", "--permPrecision", dest="permPrecision", \
      default=1, type=float,\
      help="numeric>0, default=.1, inverse of number of permutations")
  parser.add_argument("-x", "--theoPrecision", dest="theoPrecision", \
      default=0.1, type=float,\
      help="numeric>0, default=0.1, precision of theo approximation" )
  parser.add_argument("-d", "--simDelay", dest="simDelay", default=0, type=int,\
      help="sim delay effects" )
  arg_namespace = parser.parse_args()

  # parse arguments: 
  simTimes = vars(arg_namespace)['simTimes']
  delayLimit = vars(arg_namespace)['delayLimit']
  lengthSeries = vars(arg_namespace)['lengthSeries']
  resultFile = vars(arg_namespace)['resultFile']
  approxVar = vars(arg_namespace)['approxVar']
  perm_precision = vars(arg_namespace)['permPrecision']
  theo_precision = vars(arg_namespace)['theoPrecision']
  alphaValue = vars(arg_namespace)['alphaValue']
  betaValue = vars(arg_namespace)['betaValue']
  simMethod = str.split(vars(arg_namespace)['simMethod'],',')
  simDelay = vars(arg_namespace)['simDelay']
  nullDistribution = vars(arg_namespace)['nullDistribution']
  #theo_precision = 0.0001

  if vars(arg_namespace)['trendSeries'] == None:
    trendSeries = False
    x_var = approxVar
  else:
    trendSeries = True
    x_var = approxVar
    trend_threshold = int(vars(arg_namespace)['trendSeries'])

  print >>sys.stderr, "simulating...",
  start_time = time.time()

  LS_values = np.zeros(simTimes, dtype='float')
  P_theo = np.zeros(simTimes, dtype='float')
  P_perm = np.zeros(simTimes, dtype='float')
  u1 = np.zeros(simTimes, dtype='float')
  u2 = np.zeros(simTimes, dtype='float')
  v1 = np.zeros(simTimes, dtype='float')
  v2 = np.zeros(simTimes, dtype='float')
  D = np.zeros(simTimes, dtype='int')
  Al = np.zeros(simTimes, dtype='int')
  Xs = np.zeros(simTimes, dtype='int')
  Ys = np.zeros(simTimes, dtype='int')
  alpha = np.zeros(simTimes, dtype='float')
  beta = np.zeros(simTimes, dtype='float')
  P_table = lsalib.theoPvalue(Rmax=lengthSeries, Dmax=delayLimit, \
      precision=theo_precision, x_decimal=my_decimal)

  for j in range(0, simTimes):
    #generate series
    nz_x = int(np.round(lengthSeries * (1-alphaValue))) # zeros in X
    nz_y = int(np.round(lengthSeries * (1-betaValue)))  # zeros in Y
    if not trendSeries:
      nnz_x = lengthSeries - nz_x 
      nnz_y = lengthSeries - nz_y 
    else:
      nnz_x = lengthSeries - nz_x + 1
      nnz_y = lengthSeries - nz_y + 1

    #generate two arrays: OxSeries, OySeries
    if simMethod[0] == 'idn':
      OxSeries = np.random.randn(1,nnz_x)[0]
      OySeries = np.random.randn(1,nnz_y)[0]
    elif simMethod[0] == 'sin':  # to implement
      sigma = float(simMethod[1])
      x = np.array(range(0,nnz_x+simDelay), dtype='float')
      x.shape = (1,len(x))
      sin_x = np.sin(x*np.pi/6.)
      OxSeries = sin_x[:nnz_x+1]
      OySeries = sin_x[simDelay:nnz_x+simDelay+1] + np.random.randn(1,nnz_x)
    else: # to implement
      print >>sys.stderr, "simMethod not implemented"

    lengthXSeries=len(OxSeries)
    lengthYSeries=len(OySeries)
    if not trendSeries:
      xSeries = OxSeries
      ySeries = OySeries
    else:
      xSeries = lsalib.ji_calc_trend(OxSeries, lengthXSeries, trend_threshold)
      ySeries = lsalib.ji_calc_trend(OySeries, lengthYSeries, trend_threshold)

    if nullDistribution == 'yes':
      xSeries = np.array([np.random.permutation\
          (np.concatenate((xSeries,np.zeros(nz_x))))])
      ySeries = np.array([np.random.permutation\
          (np.concatenate((ySeries,np.zeros(nz_y))))])

    xSeries = np.ma.masked_invalid(xSeries)
    ySeries = np.ma.masked_invalid(ySeries)
    (u1[j], u2[j], v1[j], v2[j]) = \
        (np.mean(xSeries),np.var(xSeries),np.mean(ySeries),np.var(ySeries))
    a = 1 - (np.sum(xSeries.mask)+np.sum(xSeries==0))/float(lengthSeries)
    b = 1 - (np.sum(ySeries.mask)+np.sum(ySeries==0))/float(lengthSeries)
    (alpha[j], beta[j]) = (a,b)

    if not trendSeries:
      xSeriesData = xSeries.reshape(1,lengthSeries)
      ySeriesData = ySeries.reshape(1,lengthSeries) 
      lsa_result = lsalib.singleLSA(xSeriesData, ySeriesData, delayLimit, \
          lsalib.simpleAverage, lsalib.percentileZNormalize, True)
      LS_values[j] = np.abs( lengthSeries * lsa_result.score )
      P_perm[j] = lsalib.permuPvalue(xSeriesData, ySeriesData, delayLimit, \
          int(1/perm_precision), LS_values[j]/lengthSeries, \
          lsalib.simpleAverage, lsalib.percentileZNormalize)
    else:
      lsa_result = lsalib.singleLSA(xSeriesData, ySeriesData, delayLimit, \
          lsalib.simpleAverage, lsalib.noneNormalize, True)
      LS_values[j] = np.abs( lengthSeries * lsa_result.score )
      P_perm[j] = lsalib.permuPvalue(xSeriesData, ySeriesData, delayLimit, \
          int(1/perm_precision), LS_values[j]/lengthSeries, \
          lsalib.simpleAverage, lsalib.noneNormalize)
    P_theo[j] = lsalib.readPvalue(P_table, R=LS_values[j], N=lengthSeries, \
        x_sd=np.sqrt(approxVar), M=1, alpha=a, beta=b, x_decimal=my_decimal)
    Al_tmp = len(lsa_result.trace)
    (Xs[j], Ys[j]) = (lsa_result.trace[Al_tmp-1][0], lsa_result.trace[Al_tmp-1][1])
    D[j]=Xs[j]-Ys[j]
    Al[j]=Al_tmp
  print >>resultFile, "R\tP_theo\tP_perm\talpha\tbeta\tu1\tu2\tv1\tv2\tXs\tYs\tD\tAl"
  print >>resultFile, '\n'.join(['\t'.join( \
    [str(LS_values[i]),str(P_theo[i]),str(P_perm[i]),str(alpha[i]),str(beta[i]),\
    str(u1[i]),str(u2[i]),str(v1[i]),str(v2[i]),\
    str(Xs[i]),str(Ys[i]),str(D[i]),str(Al[i])]\
    ) for i in xrange(0,simTimes)])

  end_time = time.time()
  elapse_time = end_time - start_time
  print >>sys.stderr, "finished in %d seconds" % elapse_time

if __name__=="__main__":
  main()
  exit(0)
