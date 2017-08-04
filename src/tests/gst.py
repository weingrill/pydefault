#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Apr 21, 2017

@author: Joerg Weingrill <jweingrill@aip.de>

compare different sidereal time calculations
'''
''' Greenwich sidereal time from UT.  Based on Meeus' _Astronomical
   Algorithms_,  pp 87-88 (2nd edition).  Note that UT0 should be
   used,  the version that reflects the earth's current rotational
   state.  (UTC comes close,  but leap seconds are inserted so that
   each second can be equal in length,  rather than "stretching" each
   second as the earth slows down.  The difference is kept within one
   second,  and can be ignored for many purposes.)   '''

def green_sidereal_time( jd_ut):
    from math import floor

    jd_ut -= 2451545.0        # set relative to 2000.0 
    t_cen = jd_ut / 36525.    # convert to julian centuries 
    """ added +1 here to match result below """
    base_t = floor( jd_ut + 1)
    jd_ut -= base_t
    rval = 280.46061837 + 360.98564736629 * jd_ut + .98564736629 * base_t + \
            t_cen * t_cen * (3.87933e-4 - t_cen / 38710000.)

    # See p 84,  in Meeus:  the following should get apparent 
    # Greenwich sidereal time:
    return( rval)

def gmst(jd):
    """returns the Greenwich mean stellar time in seconds see also http://www.cv.nrao.edu/~rfisher/Ephemerides/times.html"""
    t = (jd - 2451545.0)/36525.
    #return 100.46061837 + 36000.770053608*t + 0.000387933*t**2 - (t**3/38710000.)
    return 24110.54841 + 8640184.812866*t + 0.093104*t**2 - 6.2e-6*t**3

def gmsth(jd):
    from math import floor
    
    jd0 = floor(jd- 2451545.0)
    return 6.664520 + 0.0657098244*(jd0)

jd = 2457834.5 # 11 58 46 = 4186.0

g1 = green_sidereal_time(jd)
g2 = gmst(jd)*360./86400.
from math import floor
print "green st: ", g1
print "gmst:     ", g2
print "diff", g1 - g2
print floor(1.3),floor(-1.3), floor(2.7), floor(-2.7)
#print 1./1.00273790935 = 0.99726956633