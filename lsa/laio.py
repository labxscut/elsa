import os, sys, csv, rpy2
import xml.etree.ElementTree as etree
import xml.dom.minidom

import rpy2.rlike.container as rlc
import rpy2.robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
ro.conversion.py2ri = numpy2ri
r = ro.r

def LA_Xgmml(la_table, la_size, lsaq_table, lsaq_size, title, LA_idx=4, LS_idx=3, Delay_idx=9):

  nodes = set()
  for i in xrange(1, lsaq_size+1):  
    node_x = r['''as.character''']((lsaq_table.rx(i,True)[0]))[0]
    node_y = r['''as.character''']((lsaq_table.rx(i,True)[1]))[0]
    nodes.add(node_x)
    nodes.add(node_y)
    for a in xrange(1, la_size+1):
      node_v = r['''as.character''']((la_table.rx(a,True)[0]))[0]
      node_e = r['''as.character''']((la_table.rx(a,True)[1]))[0]
      node_z = r['''as.character''']((la_table.rx(a,True)[2]))[0]
      if ((node_x==node_v)&(node_y==node_e)):
         node_m_x_y = '_'.join( ['m', node_x, node_y] )
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
    node_x = r['''as.character''']((lsa_table.rx(i,True)[0]))[0]
    node_y = r['''as.character''']((lsa_table.rx(i,True)[1]))[0]
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

    for a in xrange(1, la_size+1):
      node_v = r['''as.character''']((la_table.rx(a,True)[0]))[0]
      node_e = r['''as.character''']((la_table.rx(a,True)[1]))[0]
      node_z = r['''as.character''']((la_table.rx(a,True)[2]))[0]
      if ((node_x==node_v)&(node_y==node_e)):
        node_m_x_y = '_'.join( ['m', node_x, node_y] )
        same += 1
        if interaction == 'pdl':
           interaction1 = 'pu'
           interaction2 = 'pdl'
        if interaction == 'ndl':
           interaction1 = 'nu'
           interaction2 = 'ndl'
        if interaction == 'pdr':
           interaction1 = 'pdr'
           interaction2 = 'pu'
        if interaction == 'ndr':
           interaction1 = 'ndr'
           interaction2 = 'nu'
        if interaction == 'pu':
           interaction1 = 'pu'
           interaction2 = 'pu'
        if interaction == 'nu':
           interaction1 = 'nu'
           interaction2 = 'nu'
        if tuple(la_table.rx(a,True)[lai])[0] >= 0:
           interaction3 = 'pu'
        else:
           interaction3 = 'nu'

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
        LA_element = etree.SubElement(edge_element, 'att')
        LA_element.set('type', 'real')
        LA_element.set('name', 'LA')
        LA_element.set('value', "%.4f" % tuple(la_table.rx(a,True)[3])[0])

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
        LA_element = etree.SubElement(edge_element, 'att')
        LA_element.set('type', 'real')
        LA_element.set('name', 'LA')
        LA_element.set('value', "%.4f" % tuple(la_table.rx(a,True)[3])[0])


        edge_label = '_'.join( [node_z, interaction3, node_m_x_y] )
        edge_element = etree.SubElement(xgmml_element, 'edge')
        edge_element.set('label', edge_label )
        edge_element.set('source', node_z )
        edge_element.set('target', node_m_x_y )
        interaction_element = etree.SubElement(edge_element, 'att')
        interaction_element.set('type', 'string')
        interaction_element.set('name', 'interaction')
        interaction_element.set('value', interaction3)
        LS_element = etree.SubElement(edge_element, 'att')
        LS_element.set('type', 'real')
        LS_element.set('name', 'LS')
        LS_element.set('value', "%.4f" % tuple(lsaq_table.rx(i,True)[2])[0])
        LA_element = etree.SubElement(edge_element, 'att')
        LA_element.set('type', 'real')
        LA_element.set('name', 'LA')
        LA_element.set('value', "%.4f" % tuple(la_table.rx(a,True)[3])[0])

    if same == 0:
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
)
