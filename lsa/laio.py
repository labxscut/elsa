import os, sys, csv, rpy2
import xml.etree.ElementTree as etree
import xml.dom.minidom
import math

import rpy2.rlike.container as rlc
import rpy2.robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
ro.conversion.py2ri = numpy2ri
r = ro.r

def tolaq(la_table, la_size, title):
  la_cols = list(r['colnames'](la_table))
  laqTable = []
  for i in xrange(0,la_size+1): 
    if i < 1:
      laqTable.append( la_cols[0:] + ["tag"])
      continue
    else:
      node_x = r['''as.character''']((la_table.rx(i,True)[0]))[0]
      node_y = r['''as.character''']((la_table.rx(i,True)[1]))[0]
      node_z = r['''as.character''']((la_table.rx(i,True)[2]))[0]
      tag =  '|'.join( [title, node_x, node_y, node_z] )
      laqTable.append([node_x, node_y, node_z] + list(r['as.character'](la_table.rx(i, True)))[3:] + [tag] )
  return laqTable 

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

def writeTable( handle, table, sep='\t' ):
  """ write a table to a delimited text file
  handle - file descriptor of that opened file
  table - table to be written
  sep - spearator, '\t' for tab delimited, ',' for comma separated
  """
  #print table
  csvWriter = csv.writer(handle, delimiter=sep, escapechar='"')
  csvWriter.writerows(table)


#def LA_Xgmml3(la_table, la_size, lsaq_table, lsaq_size, title, nodeinfor_table, nodeinfor_size, LA_idx=4, LS_idx=3, Delay_idx=9):
  #laq_nodes = set()
  #laq_edges = dict()   
  #lai = LA_idx-1 #4-1
  #di = Delay_idx-1 #9-1
  #li = LS_idx-1 #3-1

  #for i in xrange(1, laq_size+1):
    #node_x = r['''as.character''']((la_table.rx(i,True)[0]))[0]
    #node_y = r['''as.character''']((la_table.rx(i,True)[1]))[0]
    #node_z = r['''as.character''']((la_table.rx(i,True)[2]))[0]

def LA_Xgmml2(la_table, la_size, lsaq_table, lsaq_size, nodeinfor_table, nodeinfor_size, nodelist_table, nodelist_size, title, LA_idx=4, LS_idx=3, Delay_idx=9):
  node_infor=dict()
  lsaq_nodes=dict()
  lsaq_edges=dict()
  missnode=set()
  nodelist_name=set()
  for i in xrange(1, nodelist_size+1):
    # s=0
    # for k in xrange(1,115):
    #   y=r['''as.character''']((nodelist_table.rx(i,True)[k]))[0]
    #   print y
    #   if y=='0':
    #      s+=1
    #   elif y=='na':
    #      s+=1
    # if s>=58: 
    n_name=r['''as.character''']((nodelist_table.rx(i,True)[0]))[0]
    nodelist_name.add(n_name)
  for i in xrange(1, nodeinfor_size+1):
    nodename=r['''as.character''']((nodeinfor_table.rx(i,True)[0]))[0]
    nodetype=r['''as.character''']((nodeinfor_table.rx(i,True)[1]))[0]
    Do=r['''as.character''']((nodeinfor_table.rx(i,True)[19]))[0]
    L=r['''as.character''']((nodeinfor_table.rx(i,True)[39]))[0]
    P=r['''as.character''']((nodeinfor_table.rx(i,True)[20]))[0]
    C=r['''as.character''']((nodeinfor_table.rx(i,True)[21]))[0]
    O=r['''as.character''']((nodeinfor_table.rx(i,True)[22]))[0]
    F=r['''as.character''']((nodeinfor_table.rx(i,True)[23]))[0]
    G=r['''as.character''']((nodeinfor_table.rx(i,True)[24]))[0]
    SP=r['''as.character''']((nodeinfor_table.rx(i,True)[25]))[0]
    SI=r['''as.character''']((nodeinfor_table.rx(i,True)[26]))[0]
    R=r['''as.character''']((nodeinfor_table.rx(i,True)[27]))[0]
    E=r['''as.character''']((nodeinfor_table.rx(i,True)[28]))[0]
    node_infor[nodename]={"nodetype":nodetype,"Domain":Do,"6Letter":L,"Phylum":P,"Class":C,"Order":O,"Family":F,"Genus":G,"Species":SP,"SilvaTag":SI,"RDPclade":R,"Ecotype":E}
     
  #construct lsaq_nodes
  #construct lsaq_edges
  #dict={key1:value1, k2:v2, k3:v3,....}
  #dict.remove(ki)
  #dict[kj]=vj
  lai = LA_idx-1 #4-1
  di = Delay_idx-1 #9-1
  li = LS_idx-1 #3-1
  for i in xrange(1, lsaq_size+1):  
    node_x = r['''as.character''']((lsaq_table.rx(i,True)[0]))[0]
    node_y = r['''as.character''']((lsaq_table.rx(i,True)[1]))[0]
    if node_x in node_infor:
         lsaq_nodes[node_x]=node_infor[node_x]
    else:
         missnode.add(node_x)
         lsaq_nodes[node_x]={"nodetype":" ", "Domain":" ", "6Letter":" "}  
    
    if node_y in node_infor:
         lsaq_nodes[node_y]=node_infor[node_y]
    else:
         missnode.add(node_y)
         lsaq_nodes[node_y]={"nodetype":" ", "Domain":" ", "6Letter":" "}  
    if tuple(lsaq_table.rx(i,True)[di])[0] > 0:
         d_code = 'dr'      #direction reta
    elif tuple(lsaq_table.rx(i,True)[di])[0] < 0:
         d_code = 'dl'      #direction lead
    else:
         d_code = 'u'
    #print(lsa_table.rx(i,True)[li])
    if tuple(lsaq_table.rx(i,True)[li])[0] >= 0:
         c_code = 'p'
    else:
         c_code = 'n'
    interaction = c_code+d_code
    LS_score = tuple(lsaq_table.rx(i,True)[2])[0]
    LS_P = tuple(lsaq_table.rx(i,True)[9])[0]
    LS_Q = tuple(lsaq_table.rx(i,True)[20])[0]
    lsaq_edges[(node_x, node_y)]=(1, {'score':LS_score, 'interaction':interaction, 'source':node_x, 'target':node_y, 'edgetype':'LS', 'Lsp':LS_P, 'Lsq':LS_Q})
  laq_nodes=lsaq_nodes
  laq_edges=lsaq_edges

  for i in xrange(1, la_size+1):
    node_x = r['''as.character''']((la_table.rx(i,True)[0]))[0]
    node_y = r['''as.character''']((la_table.rx(i,True)[1]))[0]
    node_z = r['''as.character''']((la_table.rx(i,True)[2]))[0]
    if node_z not in nodelist_name:
       pass
    else: 
       if (node_x,node_y) in laq_edges:
         x = tuple(la_table.rx(i,True)[lai])[0]
         if isinstance(x, float) and math.isnan(x):
            #LA_score = -9999
            pass
         else:
            #LA_score = tuple(la_table.rx(i,True)[3])[0]
            LA_score = x
            node_m_x_y = '_'.join( ['m', node_x, node_y] )
            laq_nodes[node_m_x_y]={"nodetype":"mid", "Domain":"na", "6Letter":"na"} 
            if node_z in node_infor:
               laq_nodes[node_z]=node_infor[node_z]
            else:
               missnode.add(node_z)
               laq_nodes[node_z]={"nodetype":" ", "Domain":" ", "6Letter":" "}  
            if tuple(la_table.rx(i,True)[lai])[0] >= 0:
               interaction_type3 = 'pu'
            else:
               interaction_type3 = 'nu'
            if laq_edges[(node_x,node_y)][1]['interaction'] == 'pdl':
               interaction_type1 = 'pu'
               interaction_type2 = 'pdr'
            elif laq_edges[(node_x,node_y)][1]['interaction'] == 'ndl':
               interaction_type1 = 'nu'
               interaction_type2 = 'ndr'
            elif laq_edges[(node_x,node_y)][1]['interaction'] == 'pdr':
               interaction_type1 = 'pdr'
               interaction_type2 = 'pu'
            elif laq_edges[(node_x,node_y)][1]['interaction'] == 'ndr':
               interaction_type1 = 'ndr'
               interaction_type2 = 'nu'
            elif laq_edges[(node_x,node_y)][1]['interaction'] == 'pu':
               interaction_type1 = 'pu'
               interaction_type2 = 'pu'
            else:
               interaction_type1 = 'nu'
               interaction_type2 = 'nu'
            LS_score = laq_edges[(node_x,node_y)][1]['score'] 
            interaction = laq_edges[(node_x,node_y)][1]['interaction']
            LS_P = laq_edges[(node_x,node_y)][1]['Lsp']	    
            LS_Q = laq_edges[(node_x,node_y)][1]['Lsq']
            LA_P = tuple(la_table.rx(i,True)[6])[0]
            LA_Q = tuple(la_table.rx(i,True)[7])[0]
            laq_edges[(node_x, node_m_x_y)]=(-1,{'Lp':'Lsp','L_p':LS_P,'Lq':'Lsq','L_q':LS_Q,'score':LS_score,'L_name':'LS','interaction':interaction_type1,'source':node_x,'target':node_m_x_y,'edgetype':'LA'})
            laq_edges[(node_y, node_m_x_y)]=(-1,{'Lp':'Lsp','L_p':LS_P,'Lq':'Lsq','L_q':LS_Q,'score':LS_score,'L_name':'LS','interaction':interaction_type2,'source':node_y,'target':node_m_x_y,'edgetype':'LA'})
            laq_edges[(node_z, node_m_x_y)]=(-1,{'Lp':'Lap','L_p':LA_P,'Lq':'Laq','L_q':LA_Q,'score':LA_score,'L_name':'LA', 'interaction':interaction_type3,'source':node_z,'target':node_m_x_y,'edgetype':'Z'})
  print "miss node_z in nodeinfor"
  print missnode
  xgmml_element=etree.Element('graph')
  xgmml_element.set('xmlns:dc', "http://purl.org/dc/elements/1.1/")
  xgmml_element.set('xmlns:xlink', "http://www.w3.org/1999/xlink")
  xgmml_element.set('xmlns:rdf', "http://www.w3.org/1999/02/22-rdf-syntax-ns#")
  xgmml_element.set('xmlns:cy', "http://www.cytoscape.org")
  xgmml_element.set('xmlns', "http://www.cs.rpi.edu/XGMML")
  xgmml_element.set('directed', "1")
  xgmml_element.set('label', "LA Network") 

  for node in laq_nodes:
    node_element = etree.SubElement(xgmml_element, 'node')
    node_element.set('id', node)
    node_element.set('name', node)
    node_element.set('label', node)
    factorName_element = etree.SubElement(node_element, 'att')
    factorName_element.set('type', 'string')
    factorName_element.set('name', 'factorName')
    factorName_element.set('value', node)
    nodetype_element = etree.SubElement(node_element, 'att')
    nodetype_element.set('type', 'string')
    nodetype_element.set('name', 'nodetype')
    nodetype_element.set('value', laq_nodes[node]["nodetype"])
    Domain_element = etree.SubElement(node_element, 'att')
    Domain_element.set('type', 'string')
    Domain_element.set('name', 'Domain')
    Domain_element.set('value', laq_nodes[node]["Domain"])
    Letter_element = etree.SubElement(node_element, 'att')
    Letter_element.set('type', 'string')
    Letter_element.set('name', '6Letter')
    Letter_element.set('value', laq_nodes[node]["6Letter"]) 


  for edge in laq_edges:     
     if  laq_edges[edge][0] < 0: 
        edge_label = '_'.join( [laq_edges[edge][1]['source'],laq_edges[edge][1]['interaction'], laq_edges[edge][1]['target']] )
        edge_element = etree.SubElement(xgmml_element, 'edge')
        edge_element.set('label', edge_label )
        edge_element.set('source', laq_edges[edge][1]['source'])
        edge_element.set('target', laq_edges[edge][1]['target'])
        interaction_element = etree.SubElement(edge_element, 'att')
        interaction_element.set('type', 'string')
        interaction_element.set('name', 'interaction')
        interaction_element.set('value', laq_edges[edge][1]['interaction'])
        score_element = etree.SubElement(edge_element, 'att')
        score_element.set('type', 'real')       
        score_element.set('name', 'score')
        score_element.set('value', "%.4f" % laq_edges[edge][1]['score'])
        edgetype_element = etree.SubElement(edge_element, 'att')
        edgetype_element.set('type', 'string')
        edgetype_element.set('name', 'edgetype')
        edgetype_element.set('value', laq_edges[edge][1]['edgetype'])
        L_element = etree.SubElement(edge_element, 'att')
        L_element.set('type', 'real')       
        L_element.set('name', laq_edges[edge][1]['L_name'])
        L_element.set('value', "%.4f" % laq_edges[edge][1]['score'])
        LP_element = etree.SubElement(edge_element, 'att')
        LP_element.set('type', 'real')       
        LP_element.set('name', laq_edges[edge][1]['Lp'])
        LP_element.set('value', "%.4f" % laq_edges[edge][1]['L_p'])
        LQ_element = etree.SubElement(edge_element, 'att')
        LQ_element.set('type', 'real')       
        LQ_element.set('name', laq_edges[edge][1]['Lq'])
        LQ_element.set('value', "%.4f" % laq_edges[edge][1]['L_q'])

     elif laq_edges[edge][0] > 0:
        edge_label = '_'.join( [laq_edges[edge][1]['source'], laq_edges[edge][1]['interaction'], laq_edges[edge][1]['target']] )
        edge_element = etree.SubElement(xgmml_element, 'edge')
        edge_element.set('label', edge_label )
        edge_element.set('source', laq_edges[edge][1]['source'])
        edge_element.set('target', laq_edges[edge][1]['target'])
        interaction_element = etree.SubElement(edge_element, 'att')
        interaction_element.set('type', 'string')
        interaction_element.set('name', 'interaction')
        interaction_element.set('value', laq_edges[edge][1]['interaction'])
        score_element = etree.SubElement(edge_element, 'att')
        score_element.set('type', 'real')       
        score_element.set('name', 'score')
        score_element.set('value', "%.4f" % laq_edges[edge][1]['score'])
        edgetype_element = etree.SubElement(edge_element, 'att')
        edgetype_element.set('type', 'string')
        edgetype_element.set('name', 'edgetype')
        edgetype_element.set('value', laq_edges[edge][1]['edgetype'])
        L_element = etree.SubElement(edge_element, 'att')
        L_element.set('type', 'real')       
        L_element.set('name', 'LS')
        L_element.set('value', "%.4f" % laq_edges[edge][1]['score'])
        LP_element = etree.SubElement(edge_element, 'att')
        LP_element.set('type', 'real')       
        LP_element.set('name', 'Lsp')
        LP_element.set('value', "%.4f" % laq_edges[edge][1]['Lsp'])
        LQ_element = etree.SubElement(edge_element, 'att')
        LQ_element.set('type', 'real')       
        LQ_element.set('name', 'Lsq')
        LQ_element.set('value', "%.4f" % laq_edges[edge][1]['Lsq'])

     else:
        pass     
  xgmml_string = etree.tostring(xgmml_element, encoding='utf-8')
  return xml.dom.minidom.parseString(xgmml_string).toprettyxml('  ')

def LA_Xgmml(la_table, la_size, lsaq_table, lsaq_size, title, LA_idx=4, LS_idx=3, Delay_idx=9):

  nodes = set()
  for i in xrange(1, lsaq_size+1):  
    node_x = r['''as.character''']((lsaq_table.rx(i,True)[0]))[0]
    node_y = r['''as.character''']((lsaq_table.rx(i,True)[1]))[0]
    node_m_x_y = '_'.join( ['m', node_x, node_y] )
    nodes.add(node_x)
    nodes.add(node_y)
    for a in xrange(1, la_size+1):
      node_v = r['''as.character''']((la_table.rx(a,True)[0]))[0]
      node_e = r['''as.character''']((la_table.rx(a,True)[1]))[0]
      node_z = r['''as.character''']((la_table.rx(a,True)[2]))[0]
      if ((node_x==node_v)&(node_y==node_e)):
         nodes.add(node_m_x_y) 
         nodes.add(node_z)
 
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

  lai = LA_idx-1 #4-1
  di = Delay_idx-1 #9-1
  li = LS_idx-1 #3-1

  for i in xrange(1, lsaq_size+1):
    node_x = r['''as.character''']((lsaq_table.rx(i,True)[0]))[0]
    node_y = r['''as.character''']((lsaq_table.rx(i,True)[1]))[0]
    node_m_x_y = '_'.join( ['m', node_x, node_y] )
    same = 0
    if tuple(lsaq_table.rx(i,True)[di])[0] > 0:
         d_code = 'dr'      #direction reta
    elif tuple(lsaq_table.rx(i,True)[di])[0] < 0:
         d_code = 'dl'      #direction lead
    else:
         d_code = 'u'
    #print(lsa_table.rx(i,True)[li])
    if tuple(lsaq_table.rx(i,True)[li])[0] >= 0:
         c_code = 'p'
    else:
         c_code = 'n'
    interaction = c_code+d_code

    if interaction == 'pdl':
         interaction1 = 'pu'
         interaction2 = 'pdl'
    elif interaction == 'ndl':
         interaction1 = 'nu'
         interaction2 = 'ndl'
    elif interaction == 'pdr':
         interaction1 = 'pdr'
         interaction2 = 'pu'
    elif interaction == 'ndr':
         interaction1 = 'ndr'
         interaction2 = 'nu'
    elif interaction == 'pu':
         interaction1 = 'pu'
         interaction2 = 'pu'
    else:
         interaction1 = 'nu'
         interaction2 = 'nu'

    for a in xrange(1, la_size+1):
      node_v = r['''as.character''']((la_table.rx(a,True)[0]))[0]
      node_e = r['''as.character''']((la_table.rx(a,True)[1]))[0]
      node_z = r['''as.character''']((la_table.rx(a,True)[2]))[0]
      if ((node_x==node_v)&(node_y==node_e)):
        same += 1
        if tuple(la_table.rx(a,True)[lai])[0] >= 0:
           interaction3 = 'pu'
        else:
           interaction3 = 'nu'

        edge_label = '_'.join( [node_z, interaction3, node_m_x_y] )
        edge_element = etree.SubElement(xgmml_element, 'edge')
        edge_element.set('label', edge_label )
        edge_element.set('source', node_z )
        edge_element.set('target', node_m_x_y )
        interaction_element = etree.SubElement(edge_element, 'att')
        interaction_element.set('type', 'string')
        interaction_element.set('name', 'interaction')
        interaction_element.set('value', interaction3)
        LA_element = etree.SubElement(edge_element, 'att')
        LA_element.set('type', 'real')       
        LA_element.set('name', 'LA')
        LA_element.set('value', "%.4f" % tuple(la_table.rx(a,True)[3])[0])


    if same > 0:
        edge_label = '_'.join( [node_x, interaction1, node_m_x_y] )
        edge_element = etree.SubElement(xgmml_element, 'edge')
        edge_element.set('label', edge_label )
        edge_element.set('source', node_x)
        edge_element.set('target', node_m_x_y )
        interaction_element = etree.SubElement(edge_element, 'att')
        interaction_element.set('type', 'string')
        interaction_element.set('name', 'interaction')
        interaction_element.set('value', interaction1)
        LS_element = etree.SubElement(edge_element, 'att')
        LS_element.set('type', 'real')
        LS_element.set('name', 'LS')
        LS_element.set('value', "%.4f" % tuple(lsaq_table.rx(i,True)[2])[0])
           
        edge_label = '_'.join( [node_m_x_y, interaction2, node_y] )
        edge_element = etree.SubElement(xgmml_element, 'edge')
        edge_element.set('label', edge_label )
        edge_element.set('source', node_m_x_y)
        edge_element.set('target', node_y )
        interaction_element = etree.SubElement(edge_element, 'att')
        interaction_element.set('type', 'string')
        interaction_element.set('name', 'interaction')
        interaction_element.set('value', interaction2)
        LS_element = etree.SubElement(edge_element, 'att')
        LS_element.set('type', 'real')
        LS_element.set('name', 'LS')
        LS_element.set('value', "%.4f" % tuple(lsaq_table.rx(i,True)[2])[0])
    

    else:
        edge_label = '_'.join( [node_x, interaction, node_y] )
        edge_element = etree.SubElement(xgmml_element, 'edge')
        edge_element.set('label', edge_label )
        edge_element.set('source', node_x )
        edge_element.set('target', node_y )
        interaction_element = etree.SubElement(edge_element, 'att')
        interaction_element.set('type', 'string')
        interaction_element.set('name', 'interaction')
        interaction_element.set('value', interaction)
        LS_element = etree.SubElement(edge_element, 'att')
        LS_element.set('type', 'real')
        LS_element.set('name', 'LS')
        LS_element.set('value', "%.4f" % tuple(lsaq_table.rx(i,True)[2])[0])
  xgmml_string = etree.tostring(xgmml_element, encoding='utf-8')
  return xml.dom.minidom.parseString(xgmml_string).toprettyxml('  ')

# filter for SIF format table
def toSif( la_table, la_size, lsaq_table, lsaq_size, nodelist_table, nodelist_size, title, LA_idx=4, LS_idx=3, Delay_idx=9):
  """ filter the lsa table for sif format which can be accepted by Cytoscape
  table - input table
  skiprows - number of rows to skip
  """
  sifTable = []
  sifTable.append( ["X", "interaction", "Y", "LS or LA", "score", "edgetype", "P", "Q"]  )
  lsaq_edges=dict()
  nodelist_name=set()
  for i in xrange(1, nodelist_size+1):
    n_name=r['''as.character''']((nodelist_table.rx(i,True)[0]))[0]
    nodelist_name.add(n_name)
     
  lai = LA_idx-1 #4-1
  di = Delay_idx-1 #9-1
  li = LS_idx-1 #3-1
  for i in xrange(1, lsaq_size+1):  
    node_x = r['''as.character''']((lsaq_table.rx(i,True)[0]))[0]
    node_y = r['''as.character''']((lsaq_table.rx(i,True)[1]))[0]
    if tuple(lsaq_table.rx(i,True)[di])[0] > 0:
         d_code = 'dr'      #direction reta
    elif tuple(lsaq_table.rx(i,True)[di])[0] < 0:
         d_code = 'dl'      #direction lead
    else:
         d_code = 'u'
    #print(lsa_table.rx(i,True)[li])
    if tuple(lsaq_table.rx(i,True)[li])[0] >= 0:
         c_code = 'p'
    else:
         c_code = 'n'
    interaction = c_code+d_code
    LS_score = tuple(lsaq_table.rx(i,True)[2])[0]
    LS_P = tuple(lsaq_table.rx(i,True)[9])[0]
    LS_Q = tuple(lsaq_table.rx(i,True)[20])[0]
    lsaq_edges[(node_x, node_y)]=(1, {'score':LS_score, 'L_name':'LS', 'interaction':interaction, 'source':node_x, 'target':node_y, 'edgetype':'LS', 'L_p':LS_P, 'L_q':LS_Q})
  laq_edges=lsaq_edges

  for i in xrange(1, la_size+1):
    node_x = r['''as.character''']((la_table.rx(i,True)[0]))[0]
    node_y = r['''as.character''']((la_table.rx(i,True)[1]))[0]
    node_z = r['''as.character''']((la_table.rx(i,True)[2]))[0]
    if node_z not in nodelist_name:
       pass
    else: 
       if (node_x,node_y) in laq_edges:
         x = tuple(la_table.rx(i,True)[lai])[0]
         if isinstance(x, float) and math.isnan(x):
            #LA_score = -9999
            pass
         else:
            #LA_score = tuple(la_table.rx(i,True)[3])[0]
            LA_score = x
            node_m_x_y = '_'.join( ['m', node_x, node_y] )
            if tuple(la_table.rx(i,True)[lai])[0] >= 0:
               interaction_type3 = 'pu'
            else:
               interaction_type3 = 'nu'
            if laq_edges[(node_x,node_y)][1]['interaction'] == 'pdl':
               interaction_type1 = 'pu'
               interaction_type2 = 'pdr'
            elif laq_edges[(node_x,node_y)][1]['interaction'] == 'ndl':
               interaction_type1 = 'nu'
               interaction_type2 = 'ndr'
            elif laq_edges[(node_x,node_y)][1]['interaction'] == 'pdr':
               interaction_type1 = 'pdr'
               interaction_type2 = 'pu'
            elif laq_edges[(node_x,node_y)][1]['interaction'] == 'ndr':
               interaction_type1 = 'ndr'
               interaction_type2 = 'nu'
            elif laq_edges[(node_x,node_y)][1]['interaction'] == 'pu':
               interaction_type1 = 'pu'
               interaction_type2 = 'pu'
            else:
               interaction_type1 = 'nu'
               interaction_type2 = 'nu'
            LS_score = laq_edges[(node_x,node_y)][1]['score'] 
            interaction = laq_edges[(node_x,node_y)][1]['interaction']
            LS_P = laq_edges[(node_x,node_y)][1]['L_p']     
            LS_Q = laq_edges[(node_x,node_y)][1]['L_q']
            LA_P = tuple(la_table.rx(i,True)[6])[0]
            LA_Q = tuple(la_table.rx(i,True)[7])[0]
            laq_edges[(node_x, node_m_x_y)]=(-1,{'Lp':'Lsp','L_p':LS_P,'Lq':'Lsq','L_q':LS_Q,'score':LS_score,'L_name':'LS','interaction':interaction_type1,'source':node_x,'target':node_m_x_y,'edgetype':'LA'})
            laq_edges[(node_y, node_m_x_y)]=(-1,{'Lp':'Lsp','L_p':LS_P,'Lq':'Lsq','L_q':LS_Q,'score':LS_score,'L_name':'LS','interaction':interaction_type2,'source':node_y,'target':node_m_x_y,'edgetype':'LA'})
            laq_edges[(node_z, node_m_x_y)]=(-1,{'Lp':'Lap','L_p':LA_P,'Lq':'Laq','L_q':LA_Q,'score':LA_score,'L_name':'LA', 'interaction':interaction_type3,'source':node_z,'target':node_m_x_y,'edgetype':'Z'})
  
  for edge in laq_edges:  
     sifTable.append( [laq_edges[edge][1]['source'],laq_edges[edge][1]['interaction'],laq_edges[edge][1]['target'],laq_edges[edge][1]['L_name'],laq_edges[edge][1]['score'],laq_edges[edge][1]['edgetype'],laq_edges[edge][1]['L_p'],laq_edges[edge][1]['L_q']]  )
   
  return sifTable

def tonewnode(la_table, la_size, lsaq_table, lsaq_size, nodeinfor_table, nodeinfor_size, nodelist_table, nodelist_size, title):
  node_cols = list(r['colnames'](nodeinfor_table))
  nodedepth = r['''as.character''']((nodeinfor_table.rx(1,True)[2]))[0]
  nodeTable = []
  for i in xrange(0,nodeinfor_size+1):
    if i<1:
      nodeTable.append( node_cols[0:] + ["tag"]  )
      #print row[li], row[di]
      continue
    else:          
      nodeID = r['''as.character''']((nodeinfor_table.rx(i,True)[0]))[0]
      nodetype = r['''as.character''']((nodeinfor_table.rx(i,True)[1]))[0]
      aa=r['''as.character''']((nodeinfor_table.rx(i,True)[3]))[0]
      ab=r['''as.character''']((nodeinfor_table.rx(i,True)[4]))[0]
      ac=r['''as.character''']((nodeinfor_table.rx(i,True)[5]))[0]
      ad=r['''as.character''']((nodeinfor_table.rx(i,True)[6]))[0]
      ae=r['''as.character''']((nodeinfor_table.rx(i,True)[7]))[0]
      af=r['''as.character''']((nodeinfor_table.rx(i,True)[8]))[0]
      ag=r['''as.character''']((nodeinfor_table.rx(i,True)[9]))[0]
      ah=r['''as.character''']((nodeinfor_table.rx(i,True)[10]))[0]
      ai=r['''as.character''']((nodeinfor_table.rx(i,True)[11]))[0]
      aj=r['''as.character''']((nodeinfor_table.rx(i,True)[12]))[0]
      ak=r['''as.character''']((nodeinfor_table.rx(i,True)[13]))[0]
      al=r['''as.character''']((nodeinfor_table.rx(i,True)[14]))[0]
      am=r['''as.character''']((nodeinfor_table.rx(i,True)[15]))[0]
      an=r['''as.character''']((nodeinfor_table.rx(i,True)[16]))[0]
      ao=r['''as.character''']((nodeinfor_table.rx(i,True)[17]))[0]
      ap=r['''as.character''']((nodeinfor_table.rx(i,True)[18]))[0]
      aq=r['''as.character''']((nodeinfor_table.rx(i,True)[19]))[0]
      ar=r['''as.character''']((nodeinfor_table.rx(i,True)[20]))[0]
      aS=r['''as.character''']((nodeinfor_table.rx(i,True)[21]))[0]
      at=r['''as.character''']((nodeinfor_table.rx(i,True)[22]))[0]
      au=r['''as.character''']((nodeinfor_table.rx(i,True)[23]))[0]
      av=r['''as.character''']((nodeinfor_table.rx(i,True)[24]))[0]
      aw=r['''as.character''']((nodeinfor_table.rx(i,True)[25]))[0]
      ax=r['''as.character''']((nodeinfor_table.rx(i,True)[26]))[0]
      ay=r['''as.character''']((nodeinfor_table.rx(i,True)[27]))[0]
      az=r['''as.character''']((nodeinfor_table.rx(i,True)[28]))[0]
      ba=r['''as.character''']((nodeinfor_table.rx(i,True)[29]))[0]
      bb=r['''as.character''']((nodeinfor_table.rx(i,True)[30]))[0]
      bc=r['''as.character''']((nodeinfor_table.rx(i,True)[31]))[0]
      bd=r['''as.character''']((nodeinfor_table.rx(i,True)[32]))[0]
      be=r['''as.character''']((nodeinfor_table.rx(i,True)[33]))[0]
      bf=r['''as.character''']((nodeinfor_table.rx(i,True)[34]))[0]
      bg=r['''as.character''']((nodeinfor_table.rx(i,True)[35]))[0]
      bh=r['''as.character''']((nodeinfor_table.rx(i,True)[36]))[0]
      bi=r['''as.character''']((nodeinfor_table.rx(i,True)[37]))[0]
      bj=r['''as.character''']((nodeinfor_table.rx(i,True)[38]))[0]
      bk=r['''as.character''']((nodeinfor_table.rx(i,True)[39]))[0]   
      nodeTable.append( [nodeID,nodetype,nodedepth,aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,aS,at,au,av,aw,ax,ay,az,ba,bb,bc,bd,be,bf,bg,bh,bi,bj,bk," "] )
  lsaq_edges=dict()
  nodelist_name=set()
  for i in xrange(1, nodelist_size+1):
    n_name=r['''as.character''']((nodelist_table.rx(i,True)[0]))[0]
    nodelist_name.add(n_name)
  for i in xrange(1, lsaq_size+1):  
    node_x = r['''as.character''']((lsaq_table.rx(i,True)[0]))[0]
    node_y = r['''as.character''']((lsaq_table.rx(i,True)[1]))[0]
    lsaq_edges[(node_x, node_y)]=(1) 
  for i in xrange(1, la_size+1):
    node_x = r['''as.character''']((la_table.rx(i,True)[0]))[0]
    node_y = r['''as.character''']((la_table.rx(i,True)[1]))[0]
    node_z = r['''as.character''']((la_table.rx(i,True)[2]))[0]
    if node_z not in nodelist_name:
       pass
    else: 
       if (node_x,node_y) in lsaq_edges:
         x = tuple(la_table.rx(i,True)[3])[0]
         if isinstance(x, float) and math.isnan(x):
            #LA_score = -9999
            pass
         else:
            #LA_score = tuple(la_table.rx(i,True)[3])[0]
            node_m_x_y = '_'.join( ['m', node_x, node_y] )            
            tag = '|'.join([title, node_x, node_y, node_z])
            nodeTable.append( [node_m_x_y, 'LA', nodedepth,'na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na','na',node_m_x_y, tag] )
  return nodeTable  
if __name__=="__main__":
  main()
  exit(0)
