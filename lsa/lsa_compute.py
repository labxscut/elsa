#!/usr/bin/env python
#lsa-compute -- computation script for LSA package to perform lsa table calculation 

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
import sys, csv, re, os, time, argparse, string, tempfile
#numeric libs
import numpy as np
import scipy as sp
#installed import
from lsa import lsalib
from lsa import lsaio
#debug import
import lsalib
import lsaio

def main():  

  # define arguments: delayLimit, fillMethod, permuNum
  parser = argparse.ArgumentParser(description="New LSA Commandline Tool")

  parser.add_argument("dataFile", metavar= "dataFile", dest="dataFile", type=argparse.FileType('r'), help="the input data file")
  parser.add_argument("resultFile", metavar= "resultFile", dest="resultFile", type=argparse.FileType('w'), help="the output result file")
  parser.add_argument("-d", "--delayLimit", dest="delayLimit", default=3, type=int, choices=range(0,6),
                    	help="specify the maximum delay possible, default: 3,\n choices: 0 to 6")
  parser.add_argument("-p", "--permuNum", dest="permuNum", default=1000, type=int, choices=[100, 200, 500, 1000, 2000],
                    	help="specify the number of permutations for p-value estimation, default: 1000,\n \
                          choices: 100, 200, 500, 1000, 2000.")
  parser.add_argument("-r", "--repNum", dest="repNum", default=1, type=int, choices=range(1,10),
                    	help="specify the number of replicates each time spot, default: 1,\n \
                          choices: 1 to 10.")
  parser.add_argument("-s", "--spotNum", dest="spotNum", default=1, type=int, choices=range(1,1000),
                    	help="specify the number of time spots, default: 1,\n \
                          choices: 1 to 1000.")
  parser.add_argument("-b", "--bootNum", dest="bootNum", default=100, type=int, choices=[100, 200, 500, 1000],
                    	help="specify the number of bootstraps for 95% confidence interval calculation, default: 100,\n \
                          choices: 100, 200, 500, 1000.")
  parser.add_argument("-t", "--transFunc", dest="transFunc", default='simple', choices=['simple', 'SD'],
                      help="specify the method to summarize replicates data, default: simple, \n \
                          choices: simple, SD \n \
                          NOTE:               \n \
                          simple: simple averaging \n \
                          SD: standard deviation weighted averaging;" )
  parser.add_argument("-f", "--fillMethod", dest="fillMethod", default='zero', choices=['zero', 'linear', 'quad', 'cubic'],
                    	help= "specify the method to fill missing, default: zero, \n \
                           choices: zero, linear, quad, cubic               \n \
                           NOTE:                                            \n \
                           zero: fill up with zeros, reliable p-Value;      \n \
                           linear: fill up with linear intropolation;       \n \
                           quad: fill up with quadratic spline;             \n \
                           cubic: fill up with cubic spline; ") 
  
  arg_namespace = parser.parse_args()

  #get arguments
  print >>sys.stderr, "lsa-compute"
  print >>sys.stderr, "copyright Li Xia, lxia@usc.edu"
  print >>sys.stderr, "learning arguments..."
  
  delayLimit = vars(arg_namespace)['delayLimit']
  fillMethod = vars(arg_namespace)['fillMethod']
  permuNum = vars(arg_namespace)['permuNum']
  dataFile = vars(arg_namespace)['dataFile']				#dataFile
  resultFile = vars(arg_namespace)['resultFile']			#resultFile
  repNum = vars(arg_namespace)['repNum']
  spotNum = vars(arg_namespace_['spotNum']
  bootNum = vars(arg_namespace)['bootNum']
  transFunc = vars(arg_namespace)['transFunc']

  #check transFunc and repNum compatibility
  if repNum <= 5 && transFunc == 'SD':
    print >>sys.stderr, "Not enough replicates for SD-weighted averaging, fall back to simple"
    transFunc = 'simple'
  if transFunc == 'simple'
    fTransform = lsalib.simpleAverage
  else:
    fTransform = lsalib.sdAverage
  
  #start timing main
  start_time = time.time()

  #input
  try:
    print >>sys.stderr, "reading raw data from dataFile..."
    rawData=np.genfromtxt( dataFile, comment='#', delimiter='\t', missing_values=['na',''], filling_values=np.nan, usecols=xrange(1,spotNum*repNum+1) )
  except:
    print >>sys.stderr, "error reading dataFile, please check the input format, spotNum and repNum \n \
                         input shall be a tab delimited txt file with '#' lines as comments and first column as factor label. \n \
                         After that, it shall have spotNum * repNum numeric cells for repNum-replicated spotNum-spotted series data. "

  factorNum = rawData.shape[0]
  cleanData=np.zeros( ( factorNum, spotNum, repNum), dtype='float' ) # (num_rows-1) x (num_cols/repNum) x (repNum)
  for i in xrange(0, factorNum):
    for j in range(0, repNum):
      cleanData[i,j] = np.array(rawData[i][arange(j,spotNum*repNum,repNum)], dtype='float')	
          #convert a section of row to a slice in rawData
          #rawData.append( [float(v) for v in row[1:]] )
          #rawData=loadtxt(dataFile, delimiter='\t', usecols=(1,), skiprows=1, converters={ 3:lambda s: float(s or 0)} )
          #rawData=loadtxt("lf.tmp", delimiter='\t', skiprows=1, converters={ 0:lambda 'na':0 } )
          #don't know howto skip a column, just skip it by explicit operation here
  print cleanData

  for i in xrange(0, factorNum):
    for j in range(0, repNum):
      cleanData[i,j] = lsalib.fillMissing( cleanData[i,j], fillMethod )
    
  #calculation
  lsaTable = lsalib.applyAnalysis( cleanData, delayLimit=delayLimit, bootNum=bootNum, permuNum=permuNum, fTransform=fTransform )

  #output

""" this part need rewrite!!!
    #normalize and transform the data
    print >>sys.stderr, "normalizing and transforming raw data..."
    normData=lsalib.dataNormalize(rawData)
    noraData=lsalib.dataNormalize(lsalib.normalTransform(normData))
    print >>sys.stderr, "computing lsa table..."
    rawResult=lsalib.newSigTestMaxDelay(noraData, delayLimit, permuNum)    
    #rawResult=lsalib.sigTestMaxDelay(noraData, delayLimit, permuNum)    
    print >>sys.stderr, "writing results..."
    resultWriter=csv.writer(outfile, delimiter="\t",quoting=csv.QUOTE_NONE)
    label=["X","Y","LS","X0","Y0","Len","Delay","pVal","corr","corPval"]
    #resultWriter.writerow(label)
    print >>sys.stderr, "labeling the result table..."
    referTable=lsaio.tryIO(dataFile, "rU")
    factorLabels=lsaio.readFirstCol(referTable)
    lsaio.closeIO(referTable)
    queryTable=[label] + rawResult
    queryTable=lsaio.labelTable(queryTable, 1, factorLabels)
    queryTable=lsaio.labelTable(queryTable, 2, factorLabels)
    resultWriter.writerows(queryTable)
"""    

  #make sure in and out works
  #print >>sys.stderr, "testing dataFile and resultFile..."
  #try:
  #    infile=open(dataFile,'rU')
  #except IOError:
  #    print >>sys.stderr, "Error: can't read file:" + dataFile
  #    exit(2)

  #try: 
  #    outfile=open(resultFile,'w')
  #except IOError:
  #    print >>sys.stderr, "Error: can't write file:" + resultFile
  #    exit(2)

  #get the rawData
#    rawData=[]
#    print >>sys.stderr, "prefiltering dataFile and filling up missing data with method " + fillMethod 
#    tmpinout = tempfile.mkstemp()[1]
#    print >>sys.stderr, "using %s as temporary file" % tmpinout
#    try:
#        lfoutfile=open(tmpinout, 'w')
#    except IOError:
#        print >>sys.stderr, "Error: can't create temporary file for data filtering"
#        exit(20)

#    csvReader=csv.reader(infile,delimiter='\t',quoting=csv.QUOTE_NONE)
#    #add safety measures for non-friendly input files
#    #check-1: a lable row and 2 data row are required or exit(1)
#    #check-2: same length for each row, or exit(2)
#    #check-3: each numerical unit is numerical or exit(3)
#    #check-4: each row 1 data point or more is required or exit(4)
#    #check-5: num_cols is dividable by repNum
#    num_cols = 0
#    num_rows = 0
#    print >>sys.stderr, "checking dataFile format..."
#    for row in csvReader:
#        num_rows += 1
#        if num_cols == 0:
#            num_cols = len(row)
#            #print >>sys.stderr, "first row cols:", len(row)
#            continue
#        if num_cols != len(row):
#            print >>sys.stderr, "Error: not all rows in the same length, quitting"
#            exit(2)
#        if len(row) <= 1:
#            print >>sys.stderr, "Error: not enough observations or formatting error, quitting"
#            exit(4)
#        for i in range(1, len(row)):
#            try:
#                #print >>sys.stderr, "converted number=", float(row[i])
#                float(row[i]) 
#            except ValueError:
#                if row[i] != 'na' and row[i] != '': 
#                    print >>sys.stderr, row[i]
#                    print >>sys.stderr, "Error: not valid data cells in input file, quitting"
#                    exit(3)
#                if not row[i]:
#                    print >>sys.stderr, "Attention: Missing data should be explicitly marked 'na' or leave blank"
#                    print >>sys.stderr, "Within row_id: ", row[0], "column_idx:", i
#                    exit(3)
#                #print >>sys.stderr, "other reasons...", row[i]
#    #print "num_rows=", num_rows      
#    if num_rows < 3:
#        print >>sys.stderr, "Error: not enough data points or format error, quitting"
#        exit(1)
#    infile.seek(0) #reset the infile in reader to original state 
#    #check-5
#    (assert (num_cols % repNum == 0 and num_cols >= repNum)) or exit(5)
#	print >>sys.stderr, "passed dataFile checking"

#    csvWriter=csv.writer(lfoutfile,delimiter='\t',quoting=csv.QUOTE_NONE)
#    lsalib.linearFillRow(csvReader, csvWriter, replicate=repNum, placeHolder="na", mode=fillMethod)
#    lfoutfile.close()

#"""	rawData=rawData.transpose() # new file format, why need this transpose?	""" # get rid of this!
    
    #(X_f, Y_f) = lsalib.apply_ftransform( (X, Y), F )						#F: transform function
    #(X_z, Y_z) = lsalib.apply_znormalize( (X_f, Y_f), N )					#N: normalize function
    #S_max = lsalib.single_lsa( (X_z, Y_z), D )								#D: delayLimit	
    #(S_low, S_up) = lsalib.bootstrapCI( (X, Y), F, N, S_max, alpha, B )	#B: bootNum, alpha: interval size
    #p = lsalib.pvalue( (X, Y), F, N, S_max, P )							#P: permuNum
    #q = lsalib.qvalue( p_set )												#p_set: set of p-values

    
  print >>sys.stderr, "finishing up..."
  outfile.close()
  print >>sys.stderr, "result is now available in " + resultFile
  end_time=time.time()
  print >>sys.stderr, "time elapsed %f seconds" % (end_time-start_time)
  print >>sys.stderr, "Thank you for using lsa-compute, byebye"

if __name__=="__main__":
  main()
  exit(0)
