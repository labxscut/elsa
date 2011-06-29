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
import sys, csv, re, os, time, optparse
#numeric libs
import numpy, scipy
from numpy import *
from scipy import *
from numpy.random import shuffle
#from scipy.io import read_array
#libs for debug use only
import pdb, string
#import hotshot
#import hprofile
#import hpstats
#assume private lsa-lib is in the parent path
#libPath = "../lsa-lib"
#sys.path.append(libPath)
#import lsalib
#import lsaio

from lsa import lsalib
from lsa import lsaio

def main():  

    parser = optparse.OptionParser("usage: %prog [options] dataFile resultFile")

    # 3 optional arguments: delayLimit, fillMethod, permuNum
    parser.add_option("-d", "--delayLimit", dest="delayLimit", default="3", choices=["0", "1", "2", "3", "4", "5", "6"],
                      help="specify the maximum delay possible, default: 3, "
                            " choices: 0 to 6")
    parser.add_option("-m", "--fillMethod", dest="fillMethod", default="zero", choices=["zero", "in", "in/ex"],
                      help="specify the method to fill missing, default: zero, "
                            " choices: zero, in, in/ex                         "
                            " NOTE:                                            "
                            " zero: fill up with zeros, reliable p-Value;        "
                            " in: fill up with linear intropolation, better alignment, reasonable p-Value;                  "
                            " in/ex: fill up with linear in/extropolation, better alignment but unreliable p-Value;") 
    parser.add_option("-p", "--permuNum", dest="permuNum", default="1000", choices=["100", "200", "500", "1000", "2000"],
                      help="specify the number of permutations for p-value estimation, default: 1000,"
                            " choices: 100, 200, 500, 1000, 2000")
    (options, args) = parser.parse_args()

    # 2 positional arguments: dataFile, resultFile
    if len(args) != 2:
        parser.error("incorrect number of arguments, use -h to get more help")

    #get the arguments
    print >>sys.stderr, "lsa-compute"
    print >>sys.stderr, "copyright Li Xia, lxia@usc.edu"
    print >>sys.stderr, "learning arguments..."
    delayLimit = int(options.delayLimit)
    fillMethod = options.fillMethod
    permuNum = int(options.permuNum)
    dataFile = args[0]
    resultFile = args[1]

    #timing
    start_time = time.time()

    #make sure in and out works
    print >>sys.stderr, "testing dataFile and resultFile..."
    try:
        infile=open(dataFile,'rU')
    except IOError:
        print >>sys.stderr, "Error: can't read file:" + dataFile
        exit(2)

    try: 
        outfile=open(resultFile,'w')
    except IOError:
        print >>sys.stderr, "Error: can't write file:" + resultFile
        exit(2)

    #get the rawData
    print >>sys.stderr, "reading raw data from dataFile..."
    global rawData
    rawData=[]
    print >>sys.stderr, "prefiltering dataFile and filling up missing data with method " + fillMethod 
    try:
        lfoutfile=open("lf.tmp",'w')
    except IOError:
        print >>sys.stderr, "Error: can't create temporary file for data filtering"
        exit(3)

    csvReader=csv.reader(infile,delimiter='\t',quoting=csv.QUOTE_NONE)
    #add safety measure for non-friendly input files
    #check-1: same length for each row
    #check-2: each numerical unit is numerical or
    #check-3: a lable row and 2 data row are required
    #check-4: each row 1 data point or more is required
    num_cols = 0
    num_rows = 0
    print >>sys.stderr, "checking dataFile format..."
    for row in csvReader:
        num_rows += 1
        if num_cols == 0:
            num_cols = len(row)
            #print >>sys.stderr, "first row cols:", len(row)
            continue
        if num_cols != len(row):
            print >>sys.stderr, "Error: not all rows in the same length, quitting"
            exit(2)
        if len(row) <= 1:
            print >>sys.stderr, "Error: not enough observations or formatting error, quitting"
            exit(2)
        for i in range(1, len(row)):
            try:
                #print >>sys.stderr, "converted number=", float(row[i])
                float(row[i]) 
            except ValueError:
                if row[i] != 'na' and row[i] != '': 
                    print >>sys.stderr, row[i]
                    print >>sys.stderr, "Error: not valid data cells in input file, quitting"
                    exit(3)
                if not row[i]:
                    print >>sys.stderr, "Attention: Missing data should be explicitly marked 'na' or leave blank"
                    print >>sys.stderr, "Within row_id: ", row[0], "column_idx:", i
                    exit(3)
                #print >>sys.stderr, "other reasons...", row[i]
    #print "num_rows=", num_rows      
    if num_rows < 3:
        print >>sys.stderr, "Error: not enough data points or format error, quitting"
        exit(3)
    infile.seek(0) #reset the infile in reader to original state 
    print >>sys.stderr, "passed dataFile checking"

    csvWriter=csv.writer(lfoutfile,delimiter='\t',quoting=csv.QUOTE_NONE)
    lsalib.linearFillRow(csvReader, csvWriter, placeHolder="na", mode=fillMethod)
    lfoutfile.close()
    try:
        lfinfile = open("lf.tmp","rU")
        csvReader=csv.reader(lfinfile, delimiter='\t', quoting=csv.QUOTE_NONE)
        i = 0
        for row in csvReader:
            i += 1
            if i == 1: continue
            rawData.append( [float(v) for v in row[1:]] )
            #rawData=loadtxt(dataFile, delimiter='\t', usecols=(1,), skiprows=1, converters={ 3:lambda s: float(s or 0)} )
            #rawData=loadtxt("lf.tmp", delimiter='\t', skiprows=1, converters={ 0:lambda 'na':0 } )
            #don't know howto skip a column, just skip it by explicit operation here
        lfinfile.close()
    except IOError:
        print >>sys.stderr, "Error: can't read temporarily created data filtering file"
        exit(3)

    try:
        #print "..."
        os.remove("lf.tmp")
    except OSError:
        print >>sys.stderr, "Error: can't delete temporarily created file lf.tmp"
        print >>sys.stderr, "do it yourself later" 

    #print rawData
    rawData=array(rawData, float)
    rawData=rawData.transpose() # new file format

    #get the labels and dates
    #print "reading factor labels and dates..."
    #infile.seek(0) #return to top
    #labelReader=csv.reader(infile,delimiter="\t",quoting=csv.QUOTE_NONE)
    #dateLabels=lsaio.readFirstRow(labelReader) # dates
    #factorLabels=lsaio.readFirstCol(labelReader) # factor labels
    #infile.close()

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
    print >>sys.stderr, "finishing up..."
    outfile.close()
    print >>sys.stderr, "result is now available in " + resultFile
    end_time=time.time()
    print >>sys.stderr, "time elapsed %f seconds" % (end_time-start_time)
    print >>sys.stderr, "Thank you for using lsa-compute, byebye"

if __name__=="__main__":
    main()
    exit(0)
