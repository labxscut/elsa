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
import optparse, sys, os, csv, re
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

if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] lsaFile")

    # 3 optional arguments: delayLimit, fillMethod, permuNum
    parser.add_option("-o", "--outputFile", dest="outputFile", default="",
                      help="specify the output file, default: on the screen")
    parser.add_option("-p", "--pValue", dest="pValue", default=0.05, type="float",
                      help="specify the highest pValue threshold to display, default: 0.05") 
    parser.add_option("-c", "--pCorr", dest="pCorr", default=1., type="float",
                      help="specify the highest Pearson Correlation threshold to display, default: 1.0")
    parser.add_option("-d", "--delayLimit", dest="delayLimit", default=3, type="int",
                      help="specify the longest delay threshhold to display, default: 3")
    parser.add_option("-r", "--referFile", dest="referFile", default="",
                      help="if specified, will display factor label instead of factor indecies")
    parser.add_option("-s", "--2sif", dest="sifFile", default="",
                      help="if specified, will also produce a SIF format file for cytoscape")
    parser.add_option("-l", "--listFactors", dest="listFactors", default="",
        help="Specify the factor of interest in a list separateb by comma: f1,f2,f3")
    (options, args) = parser.parse_args()

    # 2 positional arguments: dataFile, resultFile
    if len(args) != 1:
        parser.error("incorrect number of arguments, use -h for more help")

    #get the arguments
    print >>sys.stderr, "lsa-query version 0.1.1"
    print >>sys.stderr, "copyright Li Xia, lxia@usc.edu"
    print >>sys.stderr, "learning arguments..."
    delayLimit = options.delayLimit
    lsaFile = args[0]
    outputFile = options.outputFile
    referFile = options.referFile
    pValue = options.pValue
    pCorr = options.pCorr
    delayLimit = options.delayLimit
    sifFile = options.sifFile
    listFactors = options.listFactors
    #removeZero = options.removeZero
    print >>sys.stderr, "delayLimit=" + repr(delayLimit)
    print >>sys.stderr, "threshold pValue=" + repr(pValue)
    print >>sys.stderr, "threshold pCorr=" + repr(pCorr)

print >>sys.stderr, "testing lsaFile and outputFile..."

lsaTable = lsaio.tryIO(lsaFile, "rU")
rawTable = lsaio.readTable(lsaTable, '\t')
outTable = sys.stdout
if outputFile != "":
    outTable=lsaio.tryIO(outputFile,"w")

print >>sys.stderr, "querying the lsatable..."

queryTable=lsaio.lowPartTable(rawTable, 8, pValue)
queryTable=lsaio.lowPartTable(queryTable, 9, pCorr)
queryTable=lsaio.upPartTable(queryTable, 9, -pCorr)
queryTable=lsaio.lowPartTable(queryTable, 7, float(delayLimit))
queryTable=lsaio.upPartTable(queryTable, 7, -float(delayLimit))

print >>sys.stderr, "removing trivial case where zero vectors perfectly correlated..." 

#if removeZeor == "yes": # remove zero vector trivial cases
queryTable=lsaio.nonequalPartTable(queryTable, 3, 1.)

#if referFile != "":
#    print "labeling the result table..."
#    referTable=lsaio.tryIO(referFile, "r")
#    factorLabels=lsaio.readFirstCol(referTable)
#    lsaio.closeIO(referTable)
#    queryTable=lsaio.labelTable(queryTable, 1, factorLabels)
#    queryTable=lsaio.labelTable(queryTable, 2, factorLabels)

if listFactors != "":
    if referFile == "":
        print >>sys.stderr, "Error: can't select the list of Facotrs w/o referFile, try again..."
        exit(21)
    else:
        print >>sys.stderr, "selecting entries involve interested factors..."
        listFactors = listFactors.split(',')
        queryTable=lsaio.selectFactors(queryTable, listFactors)

print >>sys.stderr, "writing up result file..."
lsaio.writeTable(outTable, queryTable, '\t')

if sifFile != "":
    print >>sys.stderr, "filtering result for as SIF file for cytoscape..."
    sifTable=lsaio.tryIO(sifFile,"w")
    lsaio.writeTable(sifTable, lsaio.toSif(queryTable), "\t")
    lsaio.closeIO(sifTable)

print >>sys.stderr, "finishing up..."
lsaio.closeIO(lsaTable)
print >>sys.stderr, "Thank you for using lsa-query, byebye!"
lsaio.closeIO(outTable)
