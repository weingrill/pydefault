#!/usr/bin/python
# -*- coding: utf-8 -*-

from numpy import sqrt, fromfunction
from sys import exit

def dist_circle(n, xcen=None ,ycen=None):
    """
    Created on Dec 10, 2012
    Form a square array where each value is its distance to a given center.
    @author: jweingrill@aip.de

    INPUTS:
      N = either  a scalar specifying the size of the N x N square output
               array, or a 2 element vector specifying the size of the
               N x M rectangular output array.

    OPTIONAL INPUTS:
      XCEN,YCEN = Scalars designating the X,Y pixel center.  These need
               not be integers, and need not be located within the
               output image.   If not supplied then the center of the output
               image is used (XCEN = YCEN = (N-1)/2.).

    OUTPUTS:
       IM  - N by N (or M x N) floating array in which the value of each 
               pixel is equal to its distance to XCEN,YCEN
    """
    if type(n) == int:
        ny = n
        nx = n
    elif len(n) == 2:
        nx = n[0]
        ny = n[1]
    else:
        exit('ERROR - Output size parameter N must contain 1 or 2 elements')
        
    if xcen is None:
        xcen = (nx-1.)/2.
    if ycen is None:
        ycen = (ny-1.)/2.
        
    im = fromfunction(lambda x,y: sqrt((x-xcen)**2+(y-ycen)**2), (nx,ny))
    return(im)
