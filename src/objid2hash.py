#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 12, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
numspace = 2000
expspace = 1000
def objid2hash(objid):
    date = objid[0:8]
    exp = int(objid[10:14])
    number = int(objid[15:])
    #print date,exp,number
    year = int(date[:4])
    month = int(date[4:6]) - 1
    day = int(date[6:8]) - 1
    #print year,month,day
    return ((((year - 2010)*12 + month)*31 + day)*expspace + exp)*numspace + number 

def hash2objid(hashval):
    number = hashval % numspace
    #hashval -= number
    hashval /= numspace
    exp = hashval % expspace
    #hashval -= exp
    hashval /= expspace
    day = hashval % 31 + 1
    #hashval -= day
    hashval /= 31
    
    month = hashval % 12 +1
    hashval /= 12
    year = hashval + 2010
    return '%04d%02d%02dA-%04d-%04d' % (year,month,day,exp,number) 

                                
objid='20101231A-0499-0900'
print objid
hashval = objid2hash(objid)
print hashval
newobjid = hash2objid(hashval)
print newobjid
assert objid==newobjid
