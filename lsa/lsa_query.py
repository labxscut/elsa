#!/usr/bin/env python

#lsa-query.py -- script to perform query task for LSA package

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
import argparse, sys, os, csv, re
try:
  # installed 
  from lsa import lsaio
except ImportError:
  # debug
  import lsaio

def main():

  # define arguments: delayLimit, fillMethod, permuNum
  parser = argparse.ArgumentParser(description="Auxillary tool to new LSA package for querying lsa results")

  parser.add_argument("rawFile", metavar= "rawFile", type=argparse.FileType('r'), help="the raw lsa file")
  parser.add_argument("entryFile", metavar= "entryFile", type=argparse.FileType('w'), help="the query result file")

  parser.add_argument("-p", "--pValue", dest="pValue", default=0.05, type=float,
                      help="specify the pValue threshold for querying, default: 0.05") 
  parser.add_argument("-c", "--PCC", dest="PCC", default=1., type=float,
                      help="specify the highest Pearson Correlation Coefficient for querying, default: 1.0")
  parser.add_argument("-q", "--qValue", dest="qValue", default=0.05, type=float,
                      help="specify the qValue threshold for querying, default: 0.05") 
  parser.add_argument("-d", "--delayLimit", dest="delayLimit", default=3, type=int,
                      help="specify the longest delay threshhold for querying, default: 3")
  parser.add_argument("-s", "--sifFile", dest="sifFile", default="",
                      help="if specified, will also produce a SIF format file for cytoscape")
  parser.add_argument("-l", "--listFactors", dest="listFactors", default="",
                      help="query only the factors of interest in the list separated by comma: f1,f2,f3")
  arg_namespace = parser.parse_args()

  #get the arguments
  print >>sys.stderr, "lsa-query"
  print >>sys.stderr, "copyright Li Xia, lxia@usc.edu"
  print >>sys.stderr, "learning arguments..."
  
  rawFile = vars(arg_namespace)['rawFile']
  entryFile = vars(arg_namespace)['entryFile']
  pValue = vars(arg_namespace)['pValue']
  PCC = vars(arg_namespace)['PCC']
  qValue = vars(arg_namespace)['qValue']
  delayLimit = vars(arg_namespace)['delayLimit']
  sifFile = vars(arg_namespace)['sifFile']
  listFactors = vars(arg_namespace)['listFactors']

  #print >>sys.stderr, "delayLimit=" + repr(delayLimit)
  #print >>sys.stderr, "threshold pValue=" + repr(pValue)
  #print >>sys.stderr, "threshold pCorr=" + repr(pCorr)

  #print >>sys.stderr, "testing lsaFile and outputFile..."

  # rawTable is a list-of-list, where the i-th list is:
  # [ f1, f2, L-score, L_low, L_up, X_start, Y_start, length, Delay, P-value, PCC,  PCC P-val,  Qvalue  ]
  # [ 1,  2,  3,       4,     5,    6,       7,       8,      9,     10,      11,   12,         13      ]

  rawTable = lsaio.readTable(rawFile, '\t')

  print >>sys.stderr, "querying the lsatable..."

  print "pValue, qValue, PCC, delayLimit, listFactors", pValue, qValue, PCC, delayLimit, listFactors
  queryTable=rawTable
  queryTable=lsaio.lowPartTable(queryTable, 10, pValue)                       #P<=pValue
  queryTable=lsaio.lowPartTable(queryTable, 13, qValue)                       #Q<=qValue
  queryTable=lsaio.lowPartTable(queryTable, 11,  PCC)                          #|pcc|<=PCC
  queryTable=lsaio.upPartTable(queryTable,  11,  -PCC)
  queryTable=lsaio.lowPartTable(queryTable, 9, float(delayLimit))             #|d|<=D
  queryTable=lsaio.upPartTable(queryTable,  9, -float(delayLimit))

  #print >>sys.stderr, "removing trivial case where zero vectors perfectly correlated..." 

  #if removeZeor == "yes": # remove zero vector trivial cases
  #queryTable=lsaio.nonequalPartTable(queryTable, 3, 1.)

  #if referFile != "":
  #    print "labeling the result table..."
  #    referTable=lsaio.tryIO(referFile, "r")
  #    factorLabels=lsaio.readFirstCol(referTable)
  #    lsaio.closeIO(referTable)
  #    queryTable=lsaio.labelTable(queryTable, 1, factorLabels)
  #    queryTable=lsaio.labelTable(queryTable, 2, factorLabels)

  if listFactors != "":
    print >>sys.stderr, "selecting entries involve interested factors..."
    listFactors = listFactors.split(',')
    queryTable=lsaio.selectFactors(queryTable, listFactors)

  #print queryTable

  print >>sys.stderr, "writing up result file..."
  lsaio.writeTable(entryFile, queryTable, '\t')

  if sifFile != "":
    print >>sys.stderr, "filtering result as a SIF file for cytoscape..."
    lsaio.writeTable(lsaio.tryIO(sifFile,'w'), lsaio.toSif(queryTable),  '\t')

  print >>sys.stderr, "finishing up..."
  print >>sys.stderr, "Thank you for using lsa-query, byebye!"

if __name__=="__main__":
  main()
  exit(0)
