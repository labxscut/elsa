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
import optparse
import numpy
import scipy
import sys
from numpy import *
from scipy import *
from scipy.io import read_array
from numpy.random import shuffle
import csv
import re
import os
#libs for debug use only
import pdb
#import hotshot
#import hprofile
#import hpstats
#assume private lsa-lib is in the parent path
#libPath = "../lsa-lib"
#sys.path.append(libPath)
#import lsalib
from lsa import lsalib
from lsa import lsaio


if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] dataFile resultFile")

    # 3 optional arguments: delayLimit, fillMethod, permuNum
    parser.add_option("-d", "--delayLimit", dest="delayLimit", default="3", choices=["1", "2", "3", "4", "5", "6"],
                      help="specify the maximum delay possible, default: 3, "
                            " choices: 1 to 6")
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
        parser.error("incorrect number of arguments")

    #get the arguments
    print "lsa-compute version 0.1.1"
    print "copyright Li Xia, lxia@usc.edu"
    print "learning arguments..."
    delayLimit = int(options.delayLimit)
    fillMethod = options.fillMethod
    permuNum = int(options.permuNum)
    dataFile = args[0]
    resultFile = args[1]

#make sure in and out works
print "testing dataFile and resultFile..."
try:
    infile=open(dataFile,'rU')
except IOError:
    print "Error: can't read file:" + dataFile
    exit(2)

try: 
    outfile=open(resultFile,'w')
except IOError:
    print "Error: can't write file:" + resultFile
    exit(2)

#get the rawData
print "reading raw data from dataFile..."
global rawData
if fillMethod != "zero":
    print "prefiltering dataFile and filling up missing data with method linear" + fillMethod +"tropolation"
    try:
        lfoutfile=open("lf.tmp",'w')
    except IOError:
        print "Error: can't creat temporary file for data filtering"
        exit(3)
    csvReader=csv.reader(infile,delimiter='\t',quoting=csv.QUOTE_NONE)
    csvWriter=csv.writer(lfoutfile,delimiter='\t',quoting=csv.QUOTE_NONE)
    lsalib.linearFillRow(csvReader, csvWriter, placeHolder="na", mode=fillMethod)
    lfoutfile.close()
    try:
        rawData=read_array("lf.tmp", separator='\t', columns=(1,-1), lines=(1,-1), missing=0)
    except IOError:
        print "Error: can't read temporarily created data filtering file"
        exit(3)
    try:
        #print "..."
        os.remove("lf.tmp")
    except OSError:
        print "Error: can't delete temporarily created file lf.tmp"
        print "do it yourself later" 
else:
    print "prefiltering dataFile and filling up missing data with zero..."
    try:
        rawData=read_array(dataFile, separator='\t', columns=(1,-1), lines=(1,-1), missing=0)
    except IOError:
        print "Error: can't read dataFile"
rawData=rawData.transpose() # new file format

#get the labels and dates
#print "reading factor labels and dates..."
#infile.seek(0) #return to top
#labelReader=csv.reader(infile,delimiter="\t",quoting=csv.QUOTE_NONE)
#dateLabels=lsaio.readFirstRow(labelReader) # dates
#factorLabels=lsaio.readFirstCol(labelReader) # factor labels
#infile.close()

#normalize and transform the data
print "normalizing and transforming raw data..."
normData=lsalib.dataNormalize(rawData)
noraData=lsalib.dataNormalize(lsalib.normalTransform(normData))
print "computing lsa table..."
#rawResult=lsalib.newSigTestMaxDelay(noraData, delayLimit, permuNum)    
rawResult=lsalib.sigTestMaxDelay(noraData, delayLimit, permuNum)    
print "writing results..."
resultWriter=csv.writer(outfile, delimiter="\t",quoting=csv.QUOTE_NONE)
label=["X","Y","LS","X0","Y0","Len","Delay","pVal","corr","corPval"]
#resultWriter.writerow(label)
print "labeling the result table..."
referTable=lsaio.tryIO(dataFile, "rU")
factorLabels=lsaio.readFirstCol(referTable)
lsaio.closeIO(referTable)
queryTable=[label] + rawResult
queryTable=lsaio.labelTable(queryTable, 1, factorLabels)
queryTable=lsaio.labelTable(queryTable, 2, factorLabels)
resultWriter.writerows(queryTable)
print "finishing up..."
outfile.close()
print "result is now available in " + resultFile
print "Thank you for using lsa-compute, byebye"
