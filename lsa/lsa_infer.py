#!/usr/bin/env python

#lsa-infer.py -- script to perform inference for LSA package

#License: BSD

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

#public libs
import optparse, sys, csv, re, os
#numeric libs
from numpy import array
import pylab
#from scipy.io import read_array
from matplotlib.lines import Line2D
#libs for debug use only
import pdb
#import hotshot
#import hprofile
#import hpstats
#assume private lsa-lib is in the parent path
#libPath = "../lsa-lib"
#sys.path.append(libPath)
#import lsaio
from lsa import lsaio
from lsa import lsalib
from lsa import lsalibx



def main():

  __script__ = "lsa_infer"
  version_desc = lsalib.safeCmd('lsa_version')
  version_print = "%s (rev: %s) - copyright Li Charlie Xia, lcxia@scut.edu.cn" \
    % (__script__, version_desc)
  print(version_print, file=sys.stderr)

  parser = optparse.OptionParser("usage: %prog [options] referFile fLabel1 fLabel2 outFile")

  parser.add_option("-n", "--normalLize", dest="normalLize", default="",
    help="if specified, plot normalized OTU abundance, otherwise percentage")
  parser.add_option("-l", "--lag", dest="lag", default="",
    help="if specified, shift the lag of two sequences")
  (options, args) = parser.parse_args()

  if len(args) != 4:
    parser.error("incorrect number of arguments, use -h for more help")

  #get the arguments
  print("lsa-infer", file=sys.stderr)
  print("copyright Li Xia, lixia@stanford.edu", file=sys.stderr)
  print("learning arguments...", file=sys.stderr)

  normalLize = options.normalLize
  lag = options.lag
  referFile = args[0]
  fLabel1 = args[1]
  fLabel2 = args[2]
  outFile = args[3]
  
  print("testing referFile and outputFile...", file=sys.stderr)
  dataTable = lsaio.tryIO(referFile, "r")
  outFig = lsaio.tryIO(outFile, "w")
  lsaio.closeIO(outFig)
  
  print("reading raw data from referFile...", file=sys.stderr)
  global rawData
  rawData = []
  try:
    #rawData = read_array(referFile, separator='\t', columns=(1,-1), lines=(1,-1), missing=0)
    csvReader = csv.reader( open(referFile, "r"), delimiter="\t", quoting=csv.QUOTE_NONE )
    i = 0
    for row in csvReader:
        i += 1
        if i == 1: 
          timeLabels = row[1:]
        vs = []
        for v in row[1:]:
          try:
            vs.append(float(v))
          except ValueError:
            vs.append(0)
        rawData.append( vs )
        #rawData=loadtxt(dataFile, delimiter='\t', usecols=(1,), skiprows=1, converters={ 3:lambda s: float(s or 0)} )
        #rawData=loadtxt("lf.tmp", delimiter='\t', skiprows=1, converters={ 0:lambda 'na':0 } )
        #don't know howto skip a column, just skip it by explicit operation here
    del csvReader
    #csvReader.close()
  except IOError:
    print("Error: can't read dataFile", file=sys.stderr)
    exit(11)
  rawData=array(rawData, float)
  rawData = rawData.transpose() # accomodate for new file format
  #print rawData
  #print "normalizing and transforming raw data..."
  #normData = lsalib.dataNormalize(rawData)
  #noraData = lsalib.dataNormalize(lsalib.normalTransform(normData))
  #print noraData

  print("identifying querying factors...", file=sys.stderr)
  factorLabels = lsaio.readFirstCol(dataTable)
  firstIndex = -1
  try:
    firstIndex = factorLabels.index(fLabel1)
  except ValueError:
    print("first factor label"+fLabel1+" not found!!!", file=sys.stderr)
    exit(12)
  #print firstIndex
  #print rawData.shape[1]
  secondIndex = -1
  try:
    secondIndex = factorLabels.index(fLabel2)
  except ValueError:
    print("first factor label"+fLabel2+" not found!!!", file=sys.stderr)
    exit(13)
  #print secondIndex
  #print rawData.shape[1]
  #print rawData[:,0]

  print("plotting...", file=sys.stderr)
  data = rawData
  if normalLize != "":
    factor1_data = lsalib.dataNormalize(data[:,firstIndex])
    factor2_data = lsalib.dataNormalize(data[:,secondIndex])
  else:
    factor1_data = data[:,firstIndex]
    factor2_data = data[:,secondIndex]
  spots = data.shape[0]
  #print len(range(1,spots+1))
  #print len(data[:,firstIndex])
  line1 = pylab.plot(list(range(1,spots+1)), factor1_data, label=fLabel1)
  line2 = pylab.plot(list(range(1,spots+1)), factor2_data, label=fLabel2)
  pylab.legend()
  #line = Line2D(range(10), range(10), linestyle='-', marker='o')
  #pylab.legend([line1,line2], [fLabel1,fLabel2])
  #pylab.legend(line1, fLabel1)
  #pylab.legend(line2, fLabel2)
  pylab.xlabel('time spots')
  #xticks=[]
  #print "time series length=", len(timeLabels)
  #print "data length=", len(factor1_data)
  #xticks=[' '] + timeLabels
  xticks=timeLabels
  #for t in range(1,spots+1):
  #  xticks[i] = str(i)+'('+xticks[i]+')'
    #if t%2 == 1:
    #  xticks = xticks + [' ']
    #else:
    #  xticks = xticks + ["%s" % t]
  #print 'xticks=', xticks
  pylab.xticks(list(range(0,spots)), xticks, color='k')
  pylab.setp( pylab.gca().get_xticklabels(), rotation=85, size=6, verticalalignment='top', horizontalalignment='left' )
  #pylab.ax.xaxis.set_minor_locator(IndexLocator)
  pylab.title('Covarying of '+fLabel1+' and '+fLabel2)
  if normalLize != "":
    pylab.ylabel('Normalized Abundance')
  else:
    pylab.ylabel('Relative Abundance in Percentage')
  pylab.draw()
  print("saving figs...", file=sys.stderr)
  pylab.savefig( outFile, format='png', bbox_inches='tight' )
  print("finishing up...", file=sys.stderr)
  lsaio.closeIO(dataTable)
  print("Thank you for using lsa-infer, byebye!", file=sys.stderr)

if __name__=="__main__":
  main()
  exit(0)
