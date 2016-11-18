#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Nov 16, 2016

@author: Joerg Weingrill <jweingrill@aip.de>

https://www.aavso.org/cgi-bin/apass_download.pl?ra=08+13+43&dec=-05+45+00&radius=0.5&outtype=1
'''

import urllib2
import urllib

class Apass(object):
    '''
    classdocs
    '''


    def __init__(self, coordinates = None, radius = None):
        '''
        Constructor
        coordinates must be SkyCoord Object
        '''
        
        self.coordinates = coordinates
        self.radius = float(radius)
        print self.radius
    
    def loadfromfile(self, filename = None):
        pass
    
    def _coordinates2string(self):
        """
        return the coordinates as strings
        we don't need a high precision for the cluster coordinates
        """
        from astropy import units as u
        #from astropy.coordinates import SkyCoord  # @UnresolvedImport
        
        c = self.coordinates
        ra_str =  c.ra.to_string(unit=u.hourangle, # @UndefinedVariable
                                 sep='+', 
                                 precision=0)  
        dec_str =  c.dec.to_string(sep='+', 
                                   precision=0)
        return ra_str, dec_str
    
    def loadfromurl(self, outtype=1):
        '''
        coordinates as a SkyCoord Object
        radius in degrees 
        outtype must be 1 to receive a csv (0 would be html)
        '''
        data = {}
        data['ra'], data['dec'] = self._coordinates2string()
        data['radius'] = self.radius
        data['outtype'] = outtype
        url_values = urllib.urlencode(data)
        url = 'https://www.aavso.org/cgi-bin/apass_download.pl'
        #full_url = url + '?' + url_values
        full_url = url +'?ra=%(ra)s&dec=%(dec)s&radius=%(radius)s&outtype=%(outtype)s' % data
        data = urllib2.urlopen(full_url)
        results = data.read()
        print results
    
if __name__ == '__main__':
    cluster = {'NGC 2281': '6:48:17.0 41:04:42',
                'NGC 6940': '20:34:26.0 28:17:00',
                'NGC 1528': '4:15:23.0 51:12:54',
                'NGC 6819': '19:41:16.0 40:11:49',
                'NGC 2422': '7:36:35.0 -14:29:00',
                'NGC 2323': '7:02:42.0 -8:23:00',
                'NGC 6709': '18:51:18.0 10:19:06',
                'NGC 6633': '18:27:15.0 6:30:30',
                'NGC 1647': '4:45:55.0 19:06:54',
                'NGC 6791': '19:20:53.0 37:46:18',
                'NGC 2236': '6:29:39.0 6:49:48',
                'IC 4756': '18:39:00.0 5:27:00'
                 }
    
    from astropy import units as u
    from astropy.coordinates import SkyCoord  # @UnresolvedImport
    
    clustercoordinates =  SkyCoord(cluster['NGC 2281'], unit=(u.hourangle, u.deg))    
    apass = Apass(clustercoordinates, radius = 0.2)
    apass.loadfromurl()