import os, sys, csv, rpy2
import xml.etree.ElementTree as etree
import xml.dom.minidom


import rpy2.rlike.container as rlc
import rpy2.robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
ro.conversion.py2ri = numpy2ri
r = ro.r

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


def LA_Xgmml2(la_table, la_size, lsaq_table, lsaq_size, title, LA_idx=4, LS_idx=3, Delay_idx=9):
  lsaq_nodes=set()
  lsaq_edges=dict()
  #construct lsaq_nodes
  #construct lsaq_edges
  #dict={key1:value1, k2:v2, k3:v3,....}
  #dict.remove(ki)
  #dict[kj]=vj
  for i in xrange(1, lsaq_size+1):  
    lsaq_edges[(node_x, node_y)]=(i, None)
    
  laq_nodes=lsaq_nodes
  laq_edges=
  
  for i in xrnge
    edge = laq_edges.remove((node_x, 
    laq_edges[(node_x, node_m)] = (-1, {'LA':la_score, 'interaction':interaction_type1})  
    laq_edges[(node_y, node_m)] = (-1, {'LA':la_score, 'interaction':interaction_type1})  
    laq_edges[(node_z, node_m)] = (-1, {'LA':la_score, 'interaction':interaction_type1})  
    
  #create graph from laq_


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
if __name__=="__main__":
  main()
  exit(0)
