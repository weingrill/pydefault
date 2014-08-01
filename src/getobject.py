#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Aug 1, 2014

@author: jwe
'''
import sys

from datasource import DataSource
from subprocess import call

scienceobject = sys.argv[1]

wifsip = DataSource(database='wifsip', 
                   user='sro', 
                   host='pina.aip.de')

query = """SELECT path,filename,filter FROM frames, science 
WHERE frames.objid=science.objid AND object like '%s%%';""" % (scienceobject)

result = wifsip.query(query)

for r in result:
    print r[0],r[1]
    source = 'sro@pina:%s/%s.fitz' % r[0:2]
    target = sys.argv[2]
    #call(['/usr/bin/scp', source, target])
    source = '%s/%s.fitz[1]' % (sys.argv[2],r[1])
    target = '%s/%s_%s.fits' % (sys.argv[2],r[1],r[2])
    
    call(['/home/jwe/bin/imcopy', source, target])


