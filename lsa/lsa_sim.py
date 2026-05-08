#!/usr/bin/env python

import argparse, time, sys
import numpy as np

#kcut_min=100
#Rmax_min=10
my_decimal = 1    # preset x.? step size 0.1 for P_table, for paper tables

try:
  #debug import
  from . import lsalib
  #np.seterr(all='raise')
except ImportError:
  #installed import
  from lsa import lsalib

def main():
  
  __script__ = "lsa_sim"
  version_desc = lsalib.safeCmd('lsa_version')
  version_print = "%s (rev: %s) - copyright Li Charlie Xia, lcxia@scut.edu.cn" \
    % (__script__, version_desc)
  print(version_print, file=sys.stderr)

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
  parser.add_argument("-T", "--trendThresh", dest="trendThresh", default=None,\
      type=float,\
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
  parser.add_argument("-p", "--permPrecision", dest="permPrecision", \
      default=1, type=float,\
      help="numeric>0, default=.1, inverse of number of permutations")
  parser.add_argument("-x", "--theoPrecision", dest="theoPrecision", \
      default=0.0001, type=float,\
      help="numeric>0, default=0.0001, precision of theo approximation" )
  parser.add_argument("-d", "--simDelay", dest="simDelay", default=0, type=int,\
      help="sim delay effects" )
  parser.add_argument("-N", "--normMethod", dest="normMethod", default="percentileZ",\
      choices=["percentileZ", "none", "pnz", "percentile"],
      help="normalization method")
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
  trendThresh = vars(arg_namespace)['trendThresh']
  normMethod = vars(arg_namespace)['normMethod']
  #theo_precision = 0.0001

  #check normMethod
  if normMethod == 'none':
    zNormalize = lsalib.noneNormalize
  elif normMethod == 'percentile':
    zNormalize = lsalib.percentileNormalize
  elif normMethod == 'percentileZ':
    zNormalize = lsalib.percentileZNormalize
  elif normMethod == 'pnz':
    zNormalize = lsalib.noZeroNormalize
  else:
    zNormalize = lsalib.percentileZNormalize # fallback to default

  #if vars(arg_namespace)['trendThresh'] == None:
  #  trendSeries = False
  #  x_var = approxVar
  #else:
  #  trendSeries = True
  #  x_var = approxVar
  #  trend_threshold = int(vars(arg_namespace)['trendSeries'])

  print("simulating...", file=sys.stderr)
  assert trendThresh == None or trendThresh >= 0
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
  print("P table generated", file=sys.stderr)

  for j in range(0, simTimes):
    #generate series
    #nz_x = int(np.round(lengthSeries * (1-alphaValue))) # zeros in X, default=0
    #nz_y = int(np.round(lengthSeries * (1-betaValue)))  # zeros in Y, default=0
    #if trendThresh == None:
    #  nnz_x = lengthSeries - nz_x 
    #  nnz_y = lengthSeries - nz_y 
    #else:
    #  nnz_x = lengthSeries - nz_x + 1
    #  nnz_y = lengthSeries - nz_y + 1

    timespots = lengthSeries 
    #timespots as in original data, lengthSeries as in converted
    if trendThresh != None:
      timespots = timespots+1

    #generate two arrays: OxSeries, OySeries
    if simMethod[0] == 'idn':
      OxSeries = np.random.randn(1,timespots)[0]
      OySeries = np.random.randn(1,timespots)[0]
    #elif simMethod[0] == 'sin':  # to implement
    #  sigma = float(simMethod[1])
    #  x = np.array(range(0,nnz_x+simDelay), dtype='float')
    #  x.shape = (1,len(x))
    #  sin_x = np.sin(x*np.pi/6.)
    #  OxSeries = sin_x[:nnz_x+1]
    #  OySeries = sin_x[simDelay:nnz_x+simDelay+1] + np.random.randn(1,nnz_x)
    else: # to implement
      print("simMethod not implemented", file=sys.stderr)

    #if nullDistribution == 'yes':
    #  xSeries = np.array([np.random.permutation\
    #      (np.concatenate((OxSeries,np.zeros(nz_x))))])
    #  ySeries = np.array([np.random.permutation\
    #      (np.concatenate((OySeries,np.zeros(nz_y))))])

    #lengthXSeries=len(OxSeries)
    #lengthYSeries=len(OySeries)

    #if trendThresh != None:
    #  xSeries = lsalib.ji_calc_trend(OxSeries, timespots, trendThresh)
    #  ySeries = lsalib.ji_calc_trend(OySeries, timespots, trendThresh)
    #else:
    #  xSeries = OxSeries
    #  ySeries = OySeries

    xSeries = np.ma.masked_invalid(OxSeries)
    ySeries = np.ma.masked_invalid(OySeries)
    (u1[j], u2[j], v1[j], v2[j]) = \
        (np.mean(xSeries),np.var(xSeries),np.mean(ySeries),np.var(ySeries))
    a = 1 - (np.sum(xSeries.mask)+np.sum(xSeries==0))/float(timespots)
    b = 1 - (np.sum(ySeries.mask)+np.sum(ySeries==0))/float(timespots)
    (alpha[j], beta[j]) = (a,b)

    xSeriesData = xSeries.reshape(1,timespots)
    ySeriesData = ySeries.reshape(1,timespots) 

    lsa_result = lsalib.singleLSA(xSeriesData, ySeriesData, delayLimit, \
          lsalib.simpleAverage, zNormalize, trendThresh, keepTrace=True)
    Al_tmp = len(lsa_result.trace)
    if Al_tmp >= 1:
      (Xs[j], Ys[j]) = (lsa_result.trace[Al_tmp-1][0], lsa_result.trace[Al_tmp-1][1])
    else:
      (Xs[j], Ys[j]) = (-1, -1)
      #print >>sys.stderr, "impossible to align:"
      #print >>sys.stderr, "X=", \
      #    lsalib.ji_calc_trend(zNormalize(lsalib.simpleAverage(xSeriesData)),\
      #        lengthSeries, trendThresh) 
      #print >>sys.stderr, "Y=", \
      #    lsalib.ji_calc_trend(zNormalize(lsalib.simpleAverage(ySeriesData)),\
      #        lengthSeries, trendThresh) 
    #print "X=", xSeriesData
    #print "Y=", ySeriesData
    #print "Xt=", \
    #    lsalib.ji_calc_trend(zNormalize(lsalib.simpleAverage(xSeriesData)),\
    #              timespots, trendThresh)
    #print "Yt=", \
    #    lsalib.ji_calc_trend(zNormalize(lsalib.simpleAverage(ySeriesData)),\
    #              timespots, trendThresh)
    #print "length=", Al_tmp
    #print "R=", lsa_result.score*lengthSeries
    #print [ (lsa_result.trace[i][0], lsa_result.trace[i][1]) for i in range(0,Al_tmp) ]
      
    D[j]=Xs[j]-Ys[j]
    LS_values[j] = lengthSeries * np.abs(lsa_result.score)
    P_perm[j] = lsalib.permuPvalue(xSeriesData, ySeriesData, delayLimit, \
          int(1/perm_precision), np.abs(lsa_result.score), \
          lsalib.simpleAverage, zNormalize, trendThresh)
    P_theo[j] = lsalib.readPvalue(P_table, R=lengthSeries * np.abs(lsa_result.score),\
        N=lengthSeries, \
        x_sd=np.sqrt(approxVar), M=1., alpha=1., beta=1., x_decimal=my_decimal)
    Al[j]=Al_tmp

  print("R\tP_theo\tP_perm\talpha\tbeta\tu1\tu2\tv1\tv2\tXs\tYs\tD\tAl", file=resultFile)
  print('\n'.join(['\t'.join( \
    [str(LS_values[i]),str(P_theo[i]),str(P_perm[i]),str(alpha[i]),str(beta[i]),\
    str(u1[i]),str(u2[i]),str(v1[i]),str(v2[i]),\
    str(Xs[i]),str(Ys[i]),str(D[i]),str(Al[i])]\
    ) for i in range(0,simTimes)]), file=resultFile)

  end_time = time.time()
  elapse_time = end_time - start_time
  print("finished in %d seconds" % elapse_time, file=sys.stderr)

if __name__=="__main__":
  main()
  exit(0)
