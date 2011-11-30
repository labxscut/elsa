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
                              help="if specified, will generate trend series")
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
  trend_threshold = int(trendSeries)

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
      xSeries = np.zeros((1,lengthSeries), dtype='float')
      for i in xrange(0, lengthSeries):
        if OxSeries[0][i] == 0 and OxSeries[0][i+1] > 0:
          x_trend = 1
        elif OxSeries[0][i] == 0 and OxSeries[0][i+1] < 0:
          x_trend = -1
        elif OxSeries[0][i] == 0 and OxSeries[0][i+1] == 0:
          x_trend = 0
        else:
          x_trend = (OxSeries[0][i+1]-OxSeries[0][i])/np.abs(OxSeries[0][i])
        if x_trend >= trend_threshold:
          xSeries[0][i] = 1
        elif x_trend <= -trend_threshold:
          xSeries[0][i] = -1
        else:
          xSeries[0][i] = 0
      # y trend
      OySeries = np.random.randn(1,lengthSeries+1)
      ySeries = np.zeros((1,lengthSeries), dtype='float')
      for i in xrange(0, lengthSeries):
        if OySeries[0][i] == 0 and OySeries[0][i+1] > 0:
          y_trend = 1
        elif OySeries[0][i] == 0 and OySeries[0][i+1] < 0:
          y_trend = -1
        elif OySeries[0][i] == 0 and OySeries[0][i+1] == 0:
          y_trend = 0
        else:
          y_trend = (OySeries[0][i+1]-OySeries[0][i])/np.abs(OySeries[0][i])
        if y_trend >= trend_threshold:
          ySeries[0][i] = 1
        elif y_trend <= -trend_threshold:
          ySeries[0][i] = -1
        else:
          ySeries[0][i] = 0
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

