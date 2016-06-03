#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 1, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''
def airmass(height):
    """
    calulates the airmass as a function of height h in degrees
    taken fom Wikipedia, Pickering (2002)
    """
    from numpy import sin, pi,array
    
    height = array(height)
    if type(height)=='float':
        if height < 0.0: height = 0.0
        if height > 90.0: height = 90.0
    if height.ndim==1:
        height[height > 90.0] = 90.0
        height[height < 1.0] = 1.0
    
    #if min(h) < 0.0 or max(h) > 90.0:
    #    raise ValueError('h = %s must be in [0.,90.]' % h)
        
    return 1.0/sin((height + 244/(165 + 47*height**1.1))*pi/180)

def height(a):
    """
    returns the height from a given airmass a,
    since inversion of airmass function is difficult.
    """
    from scipy.optimize import newton
    from numpy import array, searchsorted
    
    # we define predefined table since Newton's method is
    # poor with bad x0 
    airtable = array([38.75, 26.64, 19.64, 15.26, 12.36, 10.33, 8.85, 7.73, 6.85, 
               6.15, 5.58, 5.1, 4.7, 4.36, 4.07, 3.81, 3.58, 3.38, 3.20,
               3.04, 2.9, 2.77, 2.65, 2.54, 2.44, 2.35, 2.27, 2.19, 2.12, 
               2.06, 1.99, 1.94, 1.88, 1.83, 1.78, 1.74, 1.7, 1.66, 1.62,
               1.59, 1.55, 1.52, 1.49, 1.46, 1.44, 1.41, 1.39, 1.37, 1.34, 
               1.32, 1.3, 1.29, 1.27, 1.25, 1.24, 1.22, 1.21, 1.19, 1.18, 
               1.17, 1.15, 1.14, 1.13, 1.12, 1.11, 1.1, 1.09, 1.09, 1.08, 
               1.07, 1.06, 1.06, 1.05, 1.05, 1.04, 1.04, 1.03, 1.03, 1.02, 
               1.02, 1.02, 1.01, 1.01, 1.01, 1.01, 1., 1., 1., 1., 1., 1.])
    
    x0 = 90. - searchsorted(airtable[::-1], a)
    f = lambda h, a: airmass(h) - a
    return newton(f, x0 = x0, args = (a,))
    
