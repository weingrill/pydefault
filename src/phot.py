#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 28, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

 objid          | character varying(19) | not null  | extended | 
 star           | integer               | not null  | plain    | 
 mag_isocor     | real                  | not null  | plain    | 
 magerr_isocor  | real                  |           | plain    | 
 mag_auto       | real                  |           | plain    | 
 magerr_auto    | real                  |           | plain    | 
 mag_aper       | real                  |           | plain    | 
 magerr_aper    | real                  |           | plain    | 
 background     | real                  |           | plain    | 
 xwin_image     | real                  |           | plain    | 
 ywin_image     | real                  |           | plain    | 
 alphawin_j2000 | double precision      |           | plain    | 
 deltawin_j2000 | double precision      |           | plain    | 
 errx2win_image | real                  |           | plain    | 
 erry2win_image | real                  |           | plain    | 
 ellipticity    | real                  |           | plain    | 
 flags          | integer               |           | plain    | 
 fwhm_image     | real                  |           | plain    | 
 class_star     | real                  |           | plain    | 
 coord          | point                 |           | plain    | 

'''
class PhotTable(list):
    '''
    classdocs
    '''


    def __init__(self, objid):
        '''
        Constructor
        '''
        from datasource import DataSource
        
        self.objid = objid
        
        self.wifsip = DataSource(database='stella', user='stella', host='pera.aip.de')
        
        keys = self.keys()
        valuearray = self.values()
        
        for value in valuearray:
            record = dict(zip(keys, value))
            self.append(record)
            
    def keys(self):
        query = """SELECT column_name, data_type, character_maximum_length
        FROM INFORMATION_SCHEMA.COLUMNS WHERE table_name = 'phot';"""
        result = self.wifsip.query(query)
        keys = [r[0] for r in result]
        return keys
    
    def values(self):
        query = """SELECT * from phot where objid = '%s'""" % self.objid
        result = self.wifsip.query(query)
        return result
       
if __name__ == '__main__':
    ph = PhotTable('20140713A-0027-0005')
    for d in ph[1]:
        print '%-20.20s' % d, ph[1][d]
    
    