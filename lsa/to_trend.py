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

  #print oSeries
  #print lengthSeries
  tSeries = np.zeros(lengthSeries-1, dtype='float')

  for i in xrange(0, lengthSeries-1):
    if oSeries[i] == 0 and oSeries[i+1] > 0:
      trend = 1
    elif oSeries[i] == 0 and oSeries[i+1] < 0:
      trend = -1
    elif oSeries[i] == 0 and oSeries[i+1] == 0:
      trend = 0
    else:
      trend = (oSeries[i+1]-oSeries[i])/np.abs(oSeries[i])


    if np.isnan(trend):
      tSeries[i] = np.nan
    elif trend > thresh:
      tSeries[i] = 1
    elif trend < -thresh:
      tSeries[i] = -1
    else:
      tSeries[i] = 0

  #print tSeries
  #exit()
  return tSeries

def main():

  # define arguments: delayLimit, fillMethod, permuNum
  parser = argparse.ArgumentParser(description="eLSA Convert Real Series to Trend Series Tool")
  parser.add_argument("dataFile", metavar= "dataFile", type=argparse.FileType('r'), help="the input data file")
  parser.add_argument("resultFile", metavar= "resultFile", type=argparse.FileType('w'), help="the output result file")
  parser.add_argument("-T", "--trendThresh", dest="trendThresh", default=0,
                              help="if specified must be a number, will generate trend series with the specified threshold")
  parser.add_argument("-S", "--spotNum", dest="spotNum", help="time spot number", default=0)
  arg_namespace = parser.parse_args()

  # parse arguments: 
  dataFile = vars(arg_namespace)['dataFile']
  resultFile = vars(arg_namespace)['resultFile']
  trend_threshold = float(vars(arg_namespace)['trendThresh'])
  spotNum = int(vars(arg_namespace)['spotNum'])

  realData=np.genfromtxt( dataFile, comments='#', delimiter='\t', missing_values=['na',''], filling_values=np.nan, usecols=range(1,spotNum+1) )
  dataFile.seek(0)  #rewind
  factorLabels=list(np.genfromtxt( dataFile, comments='#', delimiter='\t', usecols=xrange(0,1), dtype='string' ))

  start_time = time.time()

  rows, cols = realData.shape
  #print rows, cols
  #print realData
  trendData = []
  for j in range(0,rows):
      # trendify
      OSeries = realData[j,:]
     # print realData[j,:]
      TSeries = ji_calc_trend(OSeries, cols, trend_threshold)
      trendData.append(TSeries)

  #print len(trendData), len(trendData[1])
  
  print >>resultFile, '\t'.join(['#TREND'] + [str(v) for v in range(1,cols)])
  for j in range(0, rows):
    print >>resultFile, '\t'.join([factorLabels[j]] + ["na" if np.isnan(v) else "%d" % v for v in trendData[j]])

  end_time = time.time()
  elapse_time = end_time - start_time

  print >>sys.stderr, "finished in %d seconds" % elapse_time

if __name__=="__main__":
  main()
  exit(0)

