#!/usr/bin/env python

#la-query.py -- script to perform query task for LSA package

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
import argparse, sys, os, csv, re, rpy2
import numpy as np
try:
  # installed 
  from lsa import lsaio
  from lsa import la_io
except ImportError:
  # debug
  import lsaio
  import la_io
import rpy2.rlike.container as rlc
import rpy2.robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
ro.conversion.py2ri = numpy2ri
r = ro.r

#print '''setwd("%s")''' % os.environ.get('PWD')
r('''setwd("%s")''' % os.environ.get('PWD'))
#r('''library(graphics)''')
#r('''library(gplots)''')

def main():

  # define arguments: delayLimit, fillMethod, permuNum
  parser = argparse.ArgumentParser(description="Auxillary tool to new LSA package for querying la results")

  parser.add_argument("rawFile1", metavar= "rawFile1", type=argparse.FileType('rU'), help="the raw lsaq file")
  parser.add_argument("rawFile2", metavar= "rawFile2", type=argparse.FileType('rwU'), help="the raw la file")
  parser.add_argument("entryFile", metavar= "entryFile", type=argparse.FileType('w'), help="the query result file")

  parser.add_argument("-q", "--queryLine", dest="queryLine", default=None,
                      help="specify the highest pValue threshold for querying, default: None \n \
                        formatting a query: \n \
                        '[!]la$Key1[>,<,>=,<=,==,!=]V1[|,&][!]la$Key2[>,<,>=,<=,==,!=]V2[|,&]...' \n \
                        and any groupings using '(' and ')' e.g. \n \
                        '(!la$P>0.01)&(la$Q<0.01)'") 
  parser.add_argument("-x", "--xgmmlFile", dest="xgmmlFile", default="",
                      help="if specified, will also produce a XGMML format file for cytoscape")
  parser.add_argument("-s", "--sifFile", dest="sifFile", default="",
                      help="if specified, will also produce a SIF format file for backward compatibility")
  arg_namespace = parser.parse_args()

  #get the arguments
  print >>sys.stderr, "la_query ($Revision$) - copyright Li Charlie Xia, lxia@usc.edu"
  print >>sys.stderr, "learning arguments..."
  
  rawFile1 = vars(arg_namespace)['rawFile1']
  rawFile2 = vars(arg_namespace)['rawFile2']
  entryFile = vars(arg_namespace)['entryFile']
  queryLine = vars(arg_namespace)['queryLine']
  print "q=", queryLine
  xgmmlFile = vars(arg_namespace)['xgmmlFile']
  sifFile = vars(arg_namespace)['sifFile']
  analysisTitle = os.path.basename(rawFile2.name)
  rawFile1.close()
  rawFile2.ciose()
  entryFile.close()

  print >>sys.stderr, "reading the lsatable..."
  r('''lsaq <- read.delim("%s")''' % (rawFile1.name))
  r('''la <- read.delim("%s")''' % (rawFile2.name))

  try:
    print >>sys.stderr, "querying the lsatable..."
    r('''la_select <- la[%s,]''' % queryLine)
  except ValueError:
    print >>sys.stderr, "error query formatting, try again"
    quit()
  try:
    print >>sys.stderr, "writing up result file..."
    r('''write.table( la_select, file="%s", quote=FALSE, row.names=FALSE, sep='\t' )''' % entryFile.name)
  except ValueError:
    print >>sys.stderr, "no entry selected, try again"
    quit()

  
  la_size=r('''dim(la_select)''')[0]
  lsaq_size=r('''dim(lsaq)''')[0]
  #rpy2 and R interfacing debug
  #print r.lsa_select
  #print r('''dim(lsa_select)''')[0]
  #print r.lsa_select.rx(1, True)
  #print tuple(r['colnames'](r.lsa_select))
  #print tuple(r['as.character'](r.lsa_select.rx(3, True)))
  #print tuple(r.lsa_select.rx(1, True)[2])[0]
  #print r['''as.character'''](r.lsa_select.rx(1, True)[0])[0]
  #print tuple(r.lsa_select.rx(1, True)[0])

  if xgmmlFile != "":
    print >>sys.stderr, "filtering result as a XGMML file for visualization such as cytoscape..."
    print >>lsaio.tryIO(xgmmlFile,'w'), la_io.LA_Xgmml(r.la_select, la_size, r.lsaq, lsaq_size, analysisTitle)

  #if sifFile != "":
  #  print >>sys.stderr, "filtering result as a SIF file for visualization such as cytoscape..."
  #  lsaio.writeTable(lsaio.tryIO(sifFile,'w'), lsaio.toSif(r.la_select, la_size))

  print >>sys.stderr, "finishing up..."
  print >>sys.stderr, "Thank you for using lsa-query, byebye!"

if __name__=="__main__":
  main()
  exit(0)
