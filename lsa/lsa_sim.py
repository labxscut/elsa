#!/usr/bin/env python

#This simple script generates two series of n iid Normal
#input: Number of Simulations (S), Length of Sequence (L), Delay Limit (D), Output Filename (F)
#output: R

import argparse, time, sys
import numpy as np

#kcut_min=100
#Rmax_min=10
my_decimal = 3    # preset x step size for P_table
#pipi = np.pi**2 # pi^2
#pipi_inv = 1/pipi

try:
  #debug import
  import lsalib
  #np.seterr(all='raise')
except ImportError:
  #installed import
  from lsa import lsalib

def ji_calc_trend(oSeries, lengthSeries, thresh):
  #Liping Ji and Kian-Lee Tan, Bioinformatics 2005

  tSeries = np.zeros((1,lengthSeries), dtype='float')

  for i in xrange(0, lengthSeries):
    if oSeries[0][i] == 0 and oSeries[0][i+1] > 0:
      trend = 1
    elif oSeries[0][i] == 0 and oSeries[0][i+1] < 0:
      trend = -1
    elif oSeries[0][i] == 0 and oSeries[0][i+1] == 0:
      trend = 0
    else:
      trend = (oSeries[0][i+1]-oSeries[0][i])/np.abs(oSeries[0][i])

    if trend >= thresh:
      tSeries[0][i] = 1
    elif trend <= -thresh:
      tSeries[0][i] = -1
    else:
      tSeries[0][i] = 0

  #print thresh
  #print oSeries
  #print tSeries

  return tSeries

def main():

  # define arguments: delayLimit, fillMethod, permuNum
  parser = argparse.ArgumentParser(description="LSA Simulation Tool")
  parser.add_argument("resultFile", metavar= "resultFile", type=argparse.FileType('w'), help="the output result file")
  parser.add_argument("-S", "--simTimes", dest="simTimes", default=10000, type=int,
                              help="specify the times of simulation to run, default: 10000")
  parser.add_argument("-D", "--delayLimit", dest="delayLimit", default=0, type=int,
                              help="specify the maximum delay possible, default: 0")
  parser.add_argument("-L", "--lengthSeries", dest="lengthSeries", default=50, type=int,
                              help="specify the length of series to generate, default: 50")
  parser.add_argument("-T", "--trendSeries", dest="trendSeries", default=None,
                              help="if specified must be a number, will generate trend series, with the specified threshold")
  parser.add_argument("-A", "--approxVar", dest="approxVar", default=1, type=float,
                              help="numeric>0, default=1, variance of partial sum variable")

  arg_namespace = parser.parse_args()

  # parse arguments: 
  simTimes = vars(arg_namespace)['simTimes']
  delayLimit = vars(arg_namespace)['delayLimit']
  lengthSeries = vars(arg_namespace)['lengthSeries']
  resultFile = vars(arg_namespace)['resultFile']
  approxVar = vars(arg_namespace)['approxVar']
  theo_precision = 0.0001
  perm_precision = 1
  if not vars(arg_namespace)['trendSeries']:
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
  alpha = np.zeros(simTimes, dtype='float')
  beta = np.zeros(simTimes, dtype='float')
  P_table = lsalib.theoPvalue(Rmax=lengthSeries, Dmax=delayLimit, precision=theo_precision, x_decimal=my_decimal)
  for j in range(0, simTimes):
    if not trendSeries:
      xSeries = np.ma.masked_invalid(np.random.randn(1,lengthSeries))
      ySeries = np.ma.masked_invalid(np.random.randn(1,lengthSeries))
      #print xSeries.shape
      #print ySeries.shape
    else:
      # x trend
      OxSeries = np.random.randn(1,lengthSeries+1)
      xSeries = ji_calc_trend(OxSeries, lengthSeries, trend_threshold)
      # y trend
      OySeries = np.random.randn(1,lengthSeries+1)
      ySeries = ji_calc_trend(OySeries, lengthSeries, trend_threshold)
      # mask_na
      xSeries = np.ma.masked_invalid(xSeries)
      ySeries = np.ma.masked_invalid(ySeries)
      #print xSeries.shape
      #print ySeries.shape
    #singleLSA call
    (u1[j], u2[j], v1[j], v2[j]) = (np.mean(xSeries),np.var(xSeries),np.mean(ySeries),np.var(ySeries))
    LS_values[j] = np.abs( lengthSeries * lsalib.singleLSA(xSeries, ySeries, delayLimit, lsalib.simpleAverage, lsalib.noneNormalize, False).score )
    a = 1 - (np.sum(xSeries.mask)+np.sum(xSeries==0))/float(lengthSeries)
    b = 1 - (np.sum(ySeries.mask)+np.sum(ySeries==0))/float(lengthSeries)
    (alpha[j], beta[j]) = (a,b)
    #print np.sqrt(approxVar)
    P_theo[j] = lsalib.readPvalue(P_table, R=LS_values[j], N=lengthSeries, x_sd=np.sqrt(approxVar), M=1, alpha=a, beta=b, x_decimal=my_decimal)
    P_perm[j] = lsalib.permuPvalue(xSeries, ySeries, delayLimit, int(1/perm_precision), LS_values[j]/lengthSeries, lsalib.simpleAverage, lsalib.noneNormalize)
  #print LS_values

  print >>resultFile, "R\tP_theo\tP_perm\talpha\tbeta"
  print >>resultFile, '\n'.join(['\t'.join( \
    [str(LS_values[i]),str(P_theo[i]),str(P_perm[i]),str(alpha[i]),str(beta[i])]) for i in xrange(0,simTimes)])

  end_time = time.time()
  elapse_time = end_time - start_time

  print >>sys.stderr, "finished in %d seconds" % elapse_time

if __name__=="__main__":
  main()
  exit(0)

