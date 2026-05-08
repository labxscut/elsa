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
import numpy as np
try:
  # installed 
  from lsa import lsaio
  from lsa import lsalib
except ImportError:
  # debug
  from . import lsaio
  from . import lsalib
import lsa

rpy_import=False
#try:
#  import rpy2
#  import rpy2.rlike.container as rlc
#  import rpy2.robjects as ro
#  from rpy2.robjects.numpy2ri import numpy2ri
#  ro.conversion.py2ri = numpy2ri
#  r = ro.r
#  #print '''setwd("%s")''' % os.environ.get('PWD')
#  r('''setwd("%s")''' % os.environ.get('PWD'))
#  #r('''library(graphics)''')
#  #r('''library(gplots)''')
#except ImportError:
#  print >>sys.stderr, "IMPORTANT!!!: R and rpy2 are not working on this system"
#  print >>sys.stderr, "IMPORTANT!!!: This script is only workable with R and rpy2"
#  #print >>sys.stderr, "R and rpy2 is not available on this platform"
#  rpy_import=False
#  exit()

def main():

  __script__ = "lsa_query"
  version_desc = lsalib.safeCmd('lsa_version')
  version_print = "%s (rev: %s) - copyright Li Charlie Xia, lcxia@scut.edu.cn" \
    % (__script__, version_desc)
  print(version_print, file=sys.stderr)
  
  # define arguments: delayLimit, fillMethod, permuNum
  parser = argparse.ArgumentParser(description="Auxillary tool to new LSA package for querying lsa results")

  parser.add_argument("rawFile", metavar= "rawFile", type=argparse.FileType('r'), help="the raw lsa file")
  parser.add_argument("entryFile", metavar= "entryFile", type=argparse.FileType('w'), help="the query result file")

  parser.add_argument("-q", "--queryLine", dest="queryLine", default=None,
                      help="specify the highest pValue threshold for querying, default: None \n \
                        formatting a query: \n \
                        '[!]lsa$Key1[>,<,>=,<=,==,!=]V1[|,&][!]lsa$Key2[>,<,>=,<=,==,!=]V2[|,&]...' \n \
                        and any groupings using '(' and ')' e.g. \n \
                        '(!lsa$P>0.01)&(lsa$Q<0.01)', if use X to select ,must have'(lsa$X>x)&(!is.na(lsa$X))'") 
  parser.add_argument("-x", "--xgmmlFile", dest="xgmmlFile", default="",
                      help="if specified, will also produce a XGMML format file for cytoscape")
  parser.add_argument("-s", "--sifFile", dest="sifFile", default="",
                      help="if specified, will also produce a SIF format file for backward compatibility")
  arg_namespace = parser.parse_args()

  #get the arguments
  print("lsa_query ($Revision$) - copyright Li Charlie Xia, lixia@stanford.edu", file=sys.stderr)
  print("learning arguments...", file=sys.stderr)
  
  rawFile = vars(arg_namespace)['rawFile']
  entryFile = vars(arg_namespace)['entryFile']
  queryLine = vars(arg_namespace)['queryLine']
  print("q=", queryLine)
  xgmmlFile = vars(arg_namespace)['xgmmlFile']
  sifFile = vars(arg_namespace)['sifFile']
  analysisTitle = os.path.basename(rawFile.name)
  rawFile.close()
  entryFile.close()

  print("reading the lsatable...", file=sys.stderr)
  r('''lsa <- read.delim("%s")''' % rawFile.name)

  try:
    print("querying the lsatable...", file=sys.stderr)
    r('''lsa_select <- lsa[%s,]''' % queryLine)
  except ValueError:
    print("error query formatting, try again", file=sys.stderr)
    quit()
  try:
    print("writing up result file...", file=sys.stderr)
    r('''write.table( lsa_select, file="%s", quote=FALSE, row.names=FALSE, sep='\t' )''' % entryFile.name)
  except ValueError:
    print("no entry selected, try again", file=sys.stderr)
    quit()

  
  lsa_size=r('''dim(lsa_select)''')[0]
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
    print("filtering result as a XGMML file for visualization such as cytoscape...", file=sys.stderr)
    print(lsaio.toXgmml(r.lsa_select, lsa_size, analysisTitle), file=lsaio.tryIO(xgmmlFile,'w'))

  if sifFile != "":
    print("filtering result as a SIF file for visualization such as cytoscape...", file=sys.stderr)
    lsaio.writeTable(lsaio.tryIO(sifFile,'w'), lsaio.toSif(r.lsa_select, lsa_size))

  print("finishing up...", file=sys.stderr)
  print("Thank you for using lsa-query, byebye!", file=sys.stderr)

if __name__=="__main__":
  main()
  exit(0)
