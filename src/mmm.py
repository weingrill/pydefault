#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Created on Dec 6, 2012

@author: jwe
 NAME:
       MMM
 PURPOSE: 
       Estimate the sky background in a stellar contaminated field.
 EXPLANATION:  
       MMM assumes that contaminated sky pixel values overwhelmingly display 
       POSITIVE departures from the true value.  Adapted from DAOPHOT 
       routine of the same name.

 CALLING SEQUENCE:
       MMM, sky, [ skymod, sigma, skew, HIGHBAD = , READNOISE=, /DEBUG, 
                  NSKY=, /INTEGER,/SILENT]

 INPUTS:
       SKY - Array or Vector containing sky values.  This version of
               MMM does not require SKY to be sorted beforehand.  SKY
               is unaltered by this program.

 OPTIONAL OUTPUTS:
       skymod - Scalar giving estimated mode of the sky values
       SIGMA -  Scalar giving standard deviation of the peak in the sky
               histogram.  If for some reason it is impossible to derive
               skymod, then SIGMA = -1.0
       SKEW -   Scalar giving skewness of the peak in the sky histogram

               If no output variables are supplied or if /DEBUG is set
               then the values of skymod, SIGMA and SKEW will be printed.

 OPTIONAL KEYWORD INPUTS:
       HIGHBAD - scalar value of the (lowest) "bad" pixel level (e.g. cosmic 
                rays or saturated pixels) If not supplied, then there is 
                assumed to be no high bad pixels.
       MINSKY - Integer giving mininum number of sky values to be used.   MMM
                will return an error if fewer sky elements are supplied.
                Default = 20.
       MAXITER - integer giving maximum number of iterations allowed,default=50
       READNOISE - Scalar giving the read noise (or minimum noise for any 
                pixel).     Normally, MMM determines the (robust) median by 
                averaging the central 20% of the sky values.     In some cases
                where the noise is low, and pixel values are quantized a
                larger fraction may be needed.    By supplying the optional
                read noise parameter, MMM is better able to adjust the
                fraction of pixels used to determine the median.                
       /INTEGER - Set this keyword if the  input SKY vector only contains
                discrete integer values.    This keyword is only needed if the
                SKY vector is of type float or double precision, but contains 
                only discrete integer values.     (Prior to July 2004, the
                equivalent of /INTEGER was set for all data types)
       /DEBUG - If this keyword is set and non-zero, then additional 
               information is displayed at the terminal.
       /SILENT - If set, then error messages will be suppressed when MMM
                cannot compute a background.    Sigma will still be set to -1
 OPTIONAL OUTPUT KEYWORD:
      NSKY - Integer scalar giving the number of pixels actually used for the
             sky computation (after outliers have been removed).
'''
import string
from numpy import where, abs, sqrt, median, asarray, mean, std
from math import log10
from sys import exit

def mmm(sky_vector, highbad=None, debug=False,
         readnoise=0, nsky=0, discrete=False, mxiter=50, minsky=20):
# Get number of sky elements 
    nsky = len(sky_vector)            
    if nsky < minsky:
        exit('ERROR -Input vector must contain at least ' + string.strip(minsky) + ' elements')
    integer = discrete
    if not integer: 
        if type(sky_vector) is int: 
            integer = True 
# Sort SKY in ascending values
    sky = asarray(sorted(sky_vector))    
# Median value of all sky values
    skymid = median(sky)   
       
    cut1 = min([skymid - sky[0], sky[-1] - skymid]) 
    if highbad is not None: 
        cut1 = min([cut1, (highbad - skymid)])
    cut2 = skymid + cut1
    cut1 = skymid - cut1
         
# Select the pixels between Cut1 and Cut2

    good = (where(sky <= cut2) and where(sky >= cut1))[0]
    ngood = len(good) 
    if ngood == 0:
        exit('ERROR - No sky values fall within ' + cut1 + ' and ' + cut2) 
# Subtract median to improve arithmetic accuracy    
    delta = sky[good] - skymid

# Highest value accepted at upper end of vector
    maximm = max(good)
# Highest value reject at lower end of vector
    minimm = min(good) - 1               

# Compute mean and sigma (from the first pass).

# median
    skymed = median(delta) 
# mean    
    skymn = mean(delta)
# sigma
    sigma = std(delta)          
#  Add median which was subtracted off earlier
    skymn = skymn + skymid
#    If mean is less than the mode, then the contamination is slight, and the
#    mean value is what we really want.
    if (skymed < skymn):
        skymod = 3. * skymed - 2. * skymn
    else:    
        skymod = skymn

    skew = (skymn - skymod) / max([1., sigma])
    nsky = maximm - minimm 

    if debug:
        print '% MMM: Number of unrejected sky elements: ' + nsky
        print '% MMM: Mode, Sigma, Skew of sky vector:', skymod, sigma, skew   
    return([ skymod, sigma , skew, nsky])

