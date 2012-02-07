#!/usr/bin/env python

#This simple script generates two series of n iid Normal
#input: Number of Simulations (S), Length of Sequence (L), Delay Limit (D), Output Filename (F)
#output: R

import argparse, time, sys
import numpy as np

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
  arg_namespace = parser.parse_args()

  # parse arguments: 
  simTimes = vars(arg_namespace)['simTimes']
  delayLimit = vars(arg_namespace)['delayLimit']
  lengthSeries = vars(arg_namespace)['lengthSeries']
  resultFile = vars(arg_namespace)['resultFile']
  if not vars(arg_namespace)['trendSeries']:
    trendSeries = False
  else:
    trendSeries = True
    trend_threshold = int(vars(arg_namespace)['trendSeries'])

  print >>sys.stderr, "simulating...",
  start_time = time.time()

  LS_values = np.zeros(simTimes, dtype='float')
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
    LS_values[j] = lengthSeries * lsalib.singleLSA(xSeries, ySeries, delayLimit, lsalib.simpleAverage, lsalib.noneNormalize, False).score
  #print LS_values

  print >>resultFile, '\n'.join([str(v) for v in LS_values])

  end_time = time.time()
  elapse_time = end_time - start_time

  print >>sys.stderr, "finished in %d seconds" % elapse_time

if __name__=="__main__":
  main()
  exit(0)

