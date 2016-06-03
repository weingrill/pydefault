#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 20, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
def hms2dd(hms):
    """convert hours, minutes seconds to degrees"""
    import warnings
    warnings.warn('Use SkyCoords instead', DeprecationWarning)
    if type(hms) is str:
        if hms.find(':')>0: 
            hms = hms.split(':')
        else:
            hms = hms.split(' ')
    hms = [float(h) for h in hms]
    if len(hms)==2: hms.append(0.0)
    return (abs(hms[0]) + hms[1]/60. + hms[2]/3600.)*15.

def hms2hh(hms):
    """convert hours, minutes seconds to decimal hours"""
    import warnings
    warnings.warn('Use SkyCoords instead', DeprecationWarning)
    if type(hms) is str:
        if hms.find(':')>0: 
            hms = hms.split(':')
        else:
            hms = hms.split(' ')
    hms = [float(h) for h in hms]
    if len(hms)==2: hms.append(0.0)
    return (hms[0] + hms[1]/60. + hms[2]/3600.)


def dms2dd(dms):
    """convert degrees, minutes seconds to degrees"""
    from functions import sign
    import warnings
    warnings.warn('Use SkyCoords instead', DeprecationWarning)
    if type(dms) is str:
        if dms.find(':')>0: 
            dms = dms.split(':')
        else:
            dms = dms.split(' ')
    dms = [float(d) for d in dms]
    if len(dms)==2: dms.append(0.0)
    return sign(dms[0])*(abs(dms[0]) + dms[1]/60. + dms[2]/3600.)

def dd2dms(degrees):
    """convert degrees to degrees, minutes, seconds"""
    from math import trunc
    from functions import sign
    import warnings
    warnings.warn('Use SkyCoords instead', DeprecationWarning)
    adegrees = abs(degrees)
    d = trunc(adegrees)
    m = trunc((adegrees-d)*60.)
    s = ((adegrees-d)*60.-m)*60. 
    return (sign(degrees)*d,m,s)

def dd2hms(degrees):
    """convert degrees to hours, minutes, seconds"""
    from math import trunc
    import warnings
    warnings.warn('Use SkyCoords instead', DeprecationWarning)
    h = trunc(degrees/15.)
    frac_hours = degrees/15. - h
    m = trunc(frac_hours*60.)
    frac_minutes=(frac_hours)*60. - m
    s = frac_minutes*60. 
    return (h,m,s)

class Coordinates(object):
    """
    simple coordinates class.
    """
    def __init__(self, alpha, delta):
        self.alpha = alpha
        self.delta = delta
    
    def __repr__(self):
        return '(%f,%f)' % (self.alpha, self.delta)
        
    def __str__(self):
        alphadeg = dd2hms(self.alpha)
        deltadeg = dd2dms(self.delta)
        salpha = '%02d %02d %05.2f' % alphadeg
        if self.delta<0 and self.delta>-1: 
            sdelta = '-%02.2d %02d %04.1f' % deltadeg
        else:
            sdelta = '%-02.2d %02d %04.1f' % deltadeg
        return salpha + ' ' + sdelta

