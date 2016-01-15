#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 16, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import urllib2
import urllib
from HTMLParser import HTMLParser

class WebdaError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class WebdaObject(dict):
    """
    
    class to query the Webda webpage for some Open Cluster properties
    
    """
    
    def __init__(self, identifier):
        data = {}
        data['cluster'] = identifier.replace(' ','+')
        url_values = urllib.urlencode(data)
        url = 'http://www.univie.ac.at/webda/cgi-bin/ocl_page.cgi'
        full_url = url + '?' + url_values
        data = urllib2.urlopen(full_url)
        results = data.read()
        class MyHTMLParser(HTMLParser):
            data = {}
            key = None
            
            def handle_data(self, data):
                #print "data = '%s'" % data
                if self.key in ['RA','Dec']:
                    self.data[self.key] = data
                    #print '%s --> %s', (self.key,data)
                    self.key = None
                elif self.key in ['lon','lat', 'd', 'ebv', 'DM', 'age', 'FeH']:
                    try:
                        self.data[self.key] = float(data)
                    except ValueError:
                        print '%s missing' % self.key
                    #print '%s --> %f', (self.key,data)
                    self.key = None
                keytab = {'Right Ascension (2000)': 'RA',
                 'Declination (2000)': 'Dec',
                 'Galactic longitude': 'lon',
                 'Galactic latitude': 'lat',
                 'Distance [pc]': 'd',
                 'Reddening [mag]': 'ebv',
                 'Distance modulus [mag]': 'DM',
                 'Log Age': 'age',
                 'Metallicity': 'FeH'}
                if data in keytab: self.key = keytab[data]    
                
            def getdata(self):
                return self.data
                
        # instantiate the parser and feed it some HTML
        parser = MyHTMLParser()
        parser.feed(results)
        data = parser.getdata()
        #print data
        for key in data.keys():
            self[key] = data[key]
        
    def __str__(self):
        return """Right Ascension (2000): %(RA)s
                  Declination (2000):     %(Dec)s
                  Galactic longitude:     %(lon)7.3f
                  Galactic latitude:      %(lat)6.3f
                  Distance [pc]:          %(d)4.0f
                  Reddening [mag]:        %(ebv)5.3f
                  Distance modulus [mag]: %(DM)4.2f
                  Log Age:                %(age)5.3f""" % self
        
        
if __name__ == '__main__':
    wo = WebdaObject('NGC 6633')
    print wo
    #wo1 = WebdaObject('NGC 6709')
    