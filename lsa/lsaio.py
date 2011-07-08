#lsaio.py -- IO module for Local Similarity Analysis(LSA) Package

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

"""lsaio -- IO module for Local Similarity Analysis(LSA) Package

    NOTE: table here means a formatted lsa table as described in lsalib
    NOTE: since python use [0] as first element, be careful when 
    interpreting the parameter to actual list index
"""

import os, sys, csv
import xml.etree.ElementTree as etree
import xml.dom.minidom

def tryIO( file, mode ):
  """ Test the IO file before using it.
  file - string including path and filname.
  mode - input or output mode, "r" or "w", or others support by
          python open() function.
  """
  try:
    handle = open( file, mode )
  except IOError:
    print "Error: can't open " + file +" file for " + mode
    sys.exit(2)
  return handle

def closeIO( handle ):
  """ Close the IO file and release the resource.
  handle - file descriptor assigned to this IO file
  """
  try:
    handle.close()
  except IOError:
    print "Error: can't close " + handle + " , but ignored"
    
def readTable( handle, sep='\t' ):
  """ read a delimited file of a table and return a list of list 
  handle - file descriptor of that opened file
  sep - spearator, '\t' for tab delimited, ',' for comma separated
  """
  csvReader = csv.reader(handle, delimiter=sep, escapechar='"')
  table=[]
  #pdb.set_trace()
  for row in csvReader:
    table.append(row)
  return table

def readFirstLine( handle, sep='\t', startcol=1 ):
  """ read the first line of a delimited file of a table
  handle - file descriptor of that opened file
  sep - spearator, '\t' for tab delimited, ',' for comma separated
  startcol - which column to start with, default: 1, the first line
  """
  csvReader = csv.reader(handle, delimiter=sep, escapechar='"')
  line = csvReader.next() 
  return line[skipcols:]

def readFirstCol( handle, sep='\t', startrow=2 ):
  """ read the first column of a delimited file of a table
  handle - file descriptor of that opened file
  sep - spearator, '\t' for tab delimited, ',' for comma separated
  startrow - which row to start with, default: 1, the first line
  """
  csvReader = csv.reader(handle, delimiter=sep, escapechar='"')
  fcol = []
  idx = 1
  for row in csvReader:
    if idx < startrow:
      idx=idx+1
      continue
    elif row:
      fcol.append(row[0])
      #if not row[0].isspace(): #remove trailling empty lines resulting from stupid Excel conversion
  return fcol

def writeTable( handle, table, sep='\t' ):
  """ write a table to a delimited text file
  handle - file descriptor of that opened file
  table - table to be written
  sep - spearator, '\t' for tab delimited, ',' for comma separated
  """
  #print table
  csvWriter = csv.writer(handle, delimiter=sep, escapechar='"')
  csvWriter.writerows(table)

def upPartTable( table, colnum, value, skiprows=1 ):
  """ split and return a partition of table with value of 
      colnum-th colnum larger or equal than a threshold value
  table - input table
  colnum - column number of which to examine
  value - thresh hold value used to split table
  skiprows - number of rows to skip
  """
  partTable = []
  idx = 0
  for row in table:
    if idx < skiprows:
      idx =  idx + 1
      partTable.append(row)
      continue
    elif float(row[colnum-1]) >= value:
      partTable.append(row)
  return partTable

def lowPartTable( table, colnum, value, skiprows=1 ):
  """ split and return a partition of table with value of 
      colnum-th colnum smaller or equal than a threshold value
  table - input table
  colnum - column number of which to examine
  value - thresh hold value used to split table
  skiprows - number of rows to skip
  """
  partTable = []
  idx = 0
  for row in table:
    if idx < skiprows:
      idx =  idx + 1
      partTable.append(row)
      continue
    elif float(row[colnum-1]) <= value:
      partTable.append(row)
  return partTable

def equalPartTable( table, colnum, value, skiprows=1 ): 
  """ split and return a partition of table with value of 
      colnum-th colnum equal to a threshold value
  table - input table
  colnum - column number of which to examine
  value - thresh hold value used to split table
  skiprows - number of rows to skip
  """
  partTable = []
  idx = 0
  for row in table:
    if idx < skiprows:
      idx =  idx + 1
      partTable.append(row)
      continue
    elif row[colnum-1] == value:
      partTable.append(row)
  return partTable

def nonequalPartTable( table, colnum, value, skiprows=1 ):
  """ split and return a partition of table with value of 
      colnum-th colnum not equal to a specific value
  table - input table
  colnum - column number of which to examine
  value - thresh hold value used to split table
  skiprows - number of rows to skip
  """
  partTable = []
  idx = 0
  for row in table:
    if idx < skiprows:
      idx =  idx + 1
      partTable.append(row)
      continue
    elif float(row[colnum-1]) != value:
      partTable.append(row)
  return partTable

def labelTable( table, colnum, labels, skiprows=1):
  """label the colnum-th column of table with corresponding label in labels list
  table - input table
  colnum - column number of which to label, it contains the index of label in labels
  labels - a list of labels which is ordered 
  skiprows - number of rows to skip
  """
  labelTable = []
  idx = 0
  for row in table:
    if idx < skiprows:
      idx =  idx + 1
      labelTable.append(row)
      continue
    else:
      row[colnum-1]=labels[int(row[colnum-1])-1]
      labelTable.append(row)
  return labelTable

def selectFactors( table, list, skiprows=1 ):
  """ select the factor of interest to diaplay on the result
  table - input table
  list - list of interested factors
  skiprows - row to skip before start
  """
  newTable = []
  idx = 0
  for entry in table:
    if idx < skiprows:
      idx = idx + 1
      newTable.append(entry)
      continue
    elif (entry[0] in list) or (entry[1] in list):
      newTable.append(entry)
  return newTable

# filter of XGMML 
def toXgmml( table, title, skiprows=1 ):
  """ filter lsa table into xgmml format
  table - input table
  """
  
  nodes = set()
  for i in xrange(skiprows, len(table)):
    nodes.add(table[i][0])
    nodes.add(table[i][1])

  xgmml_element=etree.Element('graph')
  xgmml_element.set('xmlns:dc', "http://purl.org/dc/elements/1.1/")
  xgmml_element.set('xmlns:xlink', "http://www.w3.org/1999/xlink")
  xgmml_element.set('xmlns:rdf', "http://www.w3.org/1999/02/22-rdf-syntax-ns#")
  xgmml_element.set('xmlns:cy', "http://www.cytoscape.org")
  xgmml_element.set('xmlns', "http://www.cs.rpi.edu/XGMML")
  xgmml_element.set('directed', "1")
  xgmml_element.set('label', "LSA Network") 

  for node in nodes:
    node_element = etree.SubElement(xgmml_element, 'node')
    node_element.set('id', node)
    node_element.set('name', node)
    node_element.set('label', node)
    factorName_element = etree.SubElement(node_element, 'att')
    factorName_element.set('type', 'string')
    factorName_element.set('name', 'factorName')
    factorName_element.set('value', node)

  #[ f1, f2, LS, lowCI, upCI, Xs, Ys, Len, Delay, P, PCC,  Ppcc,  Q ]
  #[ 1,  2,  3,      4,    5,  6,  7,   8,     9,10,  11,    12, 13 ]

  di = 8
  li = 2

  for i in xrange(skiprows, len(table)):
    if table[i][di] > 0:
      rcode = 'dr'      #direction lead
    elif table[i][di] < 0:
      rcode = 'dl'      #direction retard
    else:
      rcode = 'u'
    if table[i][li] >= 0:
      ccode = 'p'
    else:
      ccode = 'n'
    interaction = ccode+rcode
    edge_label = '_'.join( [table[i][0], interaction, table[i][1]] )
    edge_element = etree.SubElement(xgmml_element, 'edge')
    edge_element.set('label', edge_label )
    edge_element.set('source', table[i][0] )
    edge_element.set('target', table[i][1] )
    interaction_element = etree.SubElement(edge_element, 'att')
    interaction_element.set('type', 'string')
    interaction_element.set('name', 'interaction')
    interaction_element.set('value', interaction)
    LS_element = etree.SubElement(edge_element, 'att')
    LS_element.set('type', 'real')
    LS_element.set('name', 'LS')
    LS_element.set('value', table[i][2])
    #P_element = etree.SubElement(edge_element, 'att')
    #P_element.set('type', 'real')
    #P_element.set('name', 'P')
    #P_element.set('value', table[i][9])
    #Q_element = etree.SubElement(edge_element, 'att')
    #Q_element.set('type','real')
    #Q_element.set('name','Q')
    #Q_element.set('value',table[i][12])

  xgmml_string = etree.tostring(xgmml_element, encoding='utf-8')
  return xml.dom.minidom.parseString(xgmml_string).toprettyxml('  ')


# filter for SIF format table
def toSif( table, LS_idx=3, Delay_idx=9, skiprows=1):
  """ filter the lsa table for sif format which can be accepted by Cytoscape
  table - input table
  skiprows - number of rows to skip
  """
  sifTable = []
  li, di = LS_idx-1, Delay_idx-1
  idx = 0
  for row in table:
    if idx < skiprows:
      idx = idx +1
      sifTable.append(["X", "interaction", "Y", row[li], row[di]])
      print row[li], row[di]
      continue
    else:
      relation = "u"
      if int(row[di]) == 0: # non delayed, undirected
        if float(row[li]) > 0:
          relation = "pu"
        elif float(row[li]) < 0:
          relation = "nu"
      elif int(row[di]) < 0: # X lead Y  
        if float(row[li]) > 0:
          relation = "pdl"
        elif float(row[li]) < 0:
          relation = "ndl"
      else: # X retard Y  
        if float(row[li]) > 0:
          relation = "pdr"
        elif float(row[li]) < 0:
          relation = "ndr"
      sifTable.append([row[0], relation, row[1], row[li], row[di]])
  return sifTable

if __name__ == "__main__":
  main()
  exit(0)
