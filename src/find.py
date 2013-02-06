#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Created on Dec 11, 2012

@author: jwe


 NAME:
    FIND
 PURPOSE:
    Find positive brightness perturbations (i.e stars) in an image 
 EXPLANATION:
    Also returns centroids and shape parameters (roundness & sharpness).
    Adapted from 1991 version of DAOPHOT, but does not allow for bad pixels
       and uses a slightly different centroid algorithm.

       Modified in March 2008 to use marginal Gaussian fits to find centroids       
 CALLING SEQUENCE:
    FIND, image, [ x, y, flux, sharp, round, hmin, fwhm, roundlim, sharplim 
        PRINT= , /SILENT, /MONITOR]

 INPUTS:
    image - 2 dimensional image array (integer or real) for which one
        wishes to identify the stars present

 OPTIONAL INPUTS:
    FIND will prompt for these parameters if not supplied

    hmin -  Threshold intensity for a point source - should generally 
        be 3 or 4 sigma above background RMS
    fwhm  - FWHM (in pixels) to be used in the convolve filter
    sharplim - 2 element vector giving low and high cutoff for the
        sharpness statistic (Default: [0.2,1.0] ).   Change this
        default only if the stars have significantly larger or 
        or smaller concentration than a Gaussian
    roundlim - 2 element vector giving low and high cutoff for the
        roundness statistic (Default: [-1.0,1.0] ).   Change this 
        default only if the stars are significantly elongated.

 OPTIONAL INPUT KEYWORDS:
       /MONITOR - Normally, FIND will display the results for each star 
                only if no output variables are supplied.   Set /MONITOR
                to always see the result of each individual star.
    /SILENT - set /SILENT keyword to suppress all output display 
    PRINT - if set and non-zero then FIND will also write its results to
        a file find.prt.   Also one can specify a different output file 
        name by setting PRINT = 'filename'.

 OPTIONAL OUTPUTS:
    x - vector containing x position of all stars identified by FIND
    y-  vector containing y position of all stars identified by FIND
    flux - vector containing flux of identified stars as determined
        by a Gaussian fit.  Fluxes are NOT converted to magnitudes.
    sharp - vector containing sharpness statistic for identified stars
    round - vector containing roundness statistic for identified stars

 NOTES:
    (1) The sharpness statistic compares the central pixel to the mean of 
       the surrounding pixels.   If this difference is greater than the 
       originally estimated height of the Gaussian or less than 0.2 the height of the
    Gaussian (for the default values of SHARPLIM) then the star will be
    rejected. 

       (2) More recent versions of FIND in DAOPHOT allow the possibility of
       ignoring bad pixels.    Unfortunately, to implement this in IDL
       would preclude the vectorization made possible with the CONVOL function
       and would run extremely slowly.

       (3) Modified in March 2008 to use marginal Gaussian distributions to 
       compute centroid.   (Formerly, find.pro determined centroids by locating
       where derivatives went to zero -- see cntrd.pro for this algorithm.   
       This was the method used in very old (~1984) versions of DAOPHOT. )   
       As discussed in more  detail in the comments to the code, the  centroid
       computation here is  the same as in IRAF DAOFIND but differs slightly 
       from the current DAOPHOT.
 PROCEDURE CALLS:
    GETOPT()
 REVISION HISTORY:
    Written W. Landsman, STX  February, 1987
    ROUND now an internal function in V3.1   W. Landsman July 1993
    Change variable name DERIV to DERIVAT    W. Landsman Feb. 1996
    Use /PRINT keyword instead of TEXTOUT    W. Landsman May  1996
    Changed loop indices to type LONG       W. Landsman Aug. 1997
       Replace DATATYPE() with size(/TNAME)   W. Landsman Nov. 2001
       Fix problem when PRINT= filename   W. Landsman   October 2002
       Fix problems with >32767 stars   D. Schlegel/W. Landsman Sep. 2004
       Fix error message when no stars found  S. Carey/W. Landsman Sep 2007
       Rewrite centroid computation to use marginal Gaussians W. Landsman 
                 Mar 2008
       Added Monitor keyword, /SILENT now suppresses all output 
                   W. Landsman    Nov 2008
       Work when threshold is negative (difference images) W. Landsman May 2010
'''

from numpy import shape, empty, where, arange, exp, fromfunction
from scipy.signal import fftconvolve
from math import sqrt
from sys import exit
import dist_circle

def gauss_kern(n, fwhm = 1.0, height = 1.0):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    n = int(n)
    c = fwhm/2.35482 # s*sqrt(2*log(2))
    c = 2 * c
    a = height
    d = dist_circle.dist_circle(n)
    g = a * exp(-( (d/c)**2 ))
    return g

def find(image,  hmin, fwhm, 
         roundlim = [-1.0, 1.0], sharplim  = [0.1, 0.9]):
    """
    find locates stars in an image
    INPUT:
        image    ... two dimensional image to search stars
        hmin     ... detection threshold
        fwhm     ... FWHM of convolution kernel in pixels
        roundlim ... roundness limits
        sharplim ... sharpness limits
    OUTPUT:
        x,y       ... coordinates of stars found that fulfill criteria
        flux      ... measured flux
        sharp     ... measured sharpness
        roundness ... measured roundness
    """
# Maximum size of convolution box in pixels
    maxbox = 13 

# Get information about the input image 
    [n_x, n_y] = shape(image)
    
    if fwhm < 0.5:
        exit('ERROR - Supplied FWHM must be at least 0.5 pixels')      
# Radius is 1.5 sigma
    radius = max([0.637*fwhm, 2.001])
    radsq = radius**2
    nhalf = min([int(radius), (maxbox-1)/2])       
# # of pixels in side of convolution box
    nbox = 2*nhalf + 1     
# Index of central pixel
    middle = nhalf          

#    lastro = n_x - nhalf
#    lastcl = n_y - nhalf
    sigsq = ( fwhm/2.35482 )**2
#Mask identifies valid pixels in convolution box    
#    mask = empty( nbox, nbox )    
#g will contain Gaussian convolution kernel
    g = gauss_kern( nbox , fwhm = fwhm)      

    row2 = (arange(nbox)-nhalf)**2

# MASK is complementary to SKIP in Stetson's Fortran
    mask = dist_circle.dist_circle([nbox, nbox]) # int(g <= radsq)     
# Value of c are now equal to distance to center
    mask = where(mask <= radsq, 1, 0)
    good = where( mask == 1)  
    pixels = sum(mask.flat)
#  Compute quantities for centroid computations that can be used for all stars

#  In fitting Gaussians to the marginal sums, pixels will arbitrarily be 
# assigned weights ranging from unity at the corners of the box to 
# NHALF^2 at the center (e.g. if NBOX = 5 or 7, the weights will be
#
#                                 1   2   3   4   3   2   1
#      1   2   3   2   1          2   4   6   8   6   4   2
#      2   4   6   4   2          3   6   9  12   9   6   3
#      3   6   9   6   3          4   8  12  16  12   8   4
#      2   4   6   4   2          3   6   9  12   9   6   3
#      1   2   3   2   1          2   4   6   8   6   4   2
#                                 1   2   3   4   3   2   1
#
# respectively).  This is done to desensitize the derived parameters to 
# possible neighboring, brighter stars.

    wt = nhalf - abs(arange(nbox)-nhalf ) + 1
    xwt = fromfunction(lambda i,j: nhalf - abs(arange(i)-nhalf ) + 1, (nbox,nbox))
    ywt = xwt.T 
    sgx = sum(g*xwt, axis = 1)
    p = sum(wt)
    sgy = sum(g*ywt, axis = 2)
    sumgx = sum(wt*sgy)
    sumgy = sum(wt*sgx)
    sumgsqy = sum(wt*sgy*sgy)
    sumgsqx = sum(wt*sgx*sgx)
    vec = nhalf - range(nbox) 
    dgdx = sgy*vec
    dgdy = sgx*vec
    sdgdxs = sum(wt*dgdx**2)
    sdgdx = sum(wt*dgdx) 
    sdgdys = sum(wt*dgdy**2)
    sdgdy = sum(wt*dgdy) 
    sgdgdx = sum(wt*sgy*dgdx)
    sgdgdy = sum(wt*sgx*dgdy)
 
# Convolution kernel now in c 
    c = g*mask               
    sumc = sum(c)
    sumcsq = sum(c**2) - sumc**2/pixels
    sumc = sumc/pixels
    c[good] = (c[good] - sumc)/sumcsq
    c1 = exp(-.5*row2/sigsq)
    sumc1 = sum(c1)/nbox
    sumc1sq = sum(c1**2) - sumc1
    c1 = (c1-sumc1)/sumc1sq

    print 'RELATIVE ERROR computed from FWHM ' + sqrt(sum(c[good]**2))

    print 'Beginning convolution of image'

# Convolve image with kernel "c"
    h = fftconvolve(image.astype(float), c, mode='same')    

    minh = h.min()
    h[:nhalf,:] = minh 
    h[n_x-nhalf:n_x,:] = minh
    h[:,:nhalf] = minh 
    h[:,n_y-nhalf:n_y] = minh

    print 'Finished convolution of image'

    mask[middle,middle] = 0    # From now on we exclude the central pixel
    pixels = pixels -1      # so the number of valid pixels is reduced by 1
    good = where(mask)      # "good" identifies position of valid pixels
    xx= (good % nbox) - middle    # x and y coordinate of valid pixels 
    yy = int(good/nbox) - middle    # relative to the center
    offset = yy*n_x + xx
#SEARCH:                 #Threshold dependent search begins here

    index = where( h >= hmin)  #Valid image pixels are greater than hmin
    nfound = len(index)
    if nfound == 0:          #Any maxima found?
        exit('ERROR - No maxima exceed input threshold of ' + hmin)

    for i in range(pixels):                             
        stars = where (h[index] >= h[index+offset[i]])
        nfound = len(stars)
# Do valid local maxima exist?
        if nfound == 0:  
            exit('ERROR - No maxima exceed input threshold of '+ hmin)
        index = index[stars]

# X index of local maxima 
    ix = index % n_x              
# Y index of local maxima
    iy = index / n_x                  
    ngood = len(index)       
    print ' local maxima located above threshold'+ngood

# NSTAR counts all stars meeting selection criteria
    nstar = 0L           
    badround = 0L 
    badsharp = 0L  
    badcntrd = 0L
# Create output X and Y arrays    
    x = empty(ngood) 
    y = empty(ngood)
# Create output flux,sharpness arrays
    flux = empty(ngood)
    sharp = empty(ngood)
    roundness = empty(ngood)

#  Loop over star positions; compute statistics

    for i in range(ngood):   
        temp = image[ix[i]-nhalf:ix[i]+nhalf,iy[i]-nhalf:iy[i]+nhalf]
# "d" is actual pixel intensity
        d = h[ix[i],iy[i]]                          

        reject = False
#  Compute Sharpness statistic

        sharp1 = (temp[middle,middle] - (sum(mask*temp))/pixels)/d
        if ( sharp1 < sharplim[0] ) or ( sharp1 > sharplim[1] ):
            badsharp = badsharp + 1
# Does not meet sharpness criteria
            reject = True

#   Compute Roundness statistic

        dx = sum( sum(temp, axis = 2)*c1)   
        dy = sum( sum(temp, axis = 1)*c1)
# Cannot compute roundness
        if (dx <= 0) or (dy <= 0):
            badround = badround + 1
            reject = True           

# Roundness statistic
        around = 2*(dx-dy) / ( dx + dy )    
# Does not meet roundness criteria
        if ( around < roundlim[0] ) or ( around > roundlim[1] ):
            badround = badround + 1
            reject = True           

#
# Centroid computation:   The centroid computation was modified in Mar 2008 and
# now differs from DAOPHOT which multiplies the correction dx by 1/(1+abs(dx)). 
# The DAOPHOT method is more robust (e.g. two different sources will not merge)
# especially in a package where the centroid will be subsequently be 
# redetermined using PSF fitting.   However, it is less accurate, and introduces
# biases in the centroid histogram.   The change here is the same made in the 
# IRAF DAOFIND routine (see 
# http://iraf.net/article.php?story=7211&query=daofind )
#    

        sd = sum(temp*ywt, axis=2)

        sumgd = sum(wt*sgy*sd)
        sumd = sum(wt*sd)
        sddgdx = sum(wt*sd*dgdx)

        hx = (sumgd - sumgx*sumd/p) / (sumgsqy - sumgx**2/p)

# HX is the height of the best-fitting marginal Gaussian.   If this is not
# positive then the centroid does not make sense 

        if (hx <= 0):
            badcntrd = badcntrd + 1
            reject = True


        skylvl = (sumd - hx*sumgx)/p
        dx = (sgdgdx - (sddgdx-sdgdx*(hx*sumgx + skylvl*p)))/(hx*sdgdxs/sigsq)
        if abs(dx) >= nhalf: 
            badcntrd = badcntrd + 1
            reject = True
# X centroid in original array
        xcen = ix[i] + dx    

# Find Y centroid                 

        sd = sum(temp*xwt, axis=1)
 
        sumgd = sum(wt*sgx*sd)
        sumd = sum(wt*sd)

        sddgdy = sum(wt*sd*dgdy)

        hy = (sumgd - sumgy*sumd/p) / (sumgsqx - sumgy**2/p)

        if (hy <= 0):
            badcntrd = badcntrd + 1
            reject = True

        skylvl = (sumd - hy*sumgy)/p
        dy = (sgdgdy - (sddgdy-sdgdy*(hy*sumgy + skylvl*p)))/(hy*sdgdys/sigsq)
        if abs(dy) >= nhalf: 
            badcntrd = badcntrd + 1
            reject = True
      
#Y centroid in original array
        ycen = iy[i] +dy    
 

#  This star has met all selection criteria.  Print out and save results

        if ~reject:
            x[nstar] = xcen 
            y[nstar] = ycen
            flux[nstar] = d 
            sharp[nstar] = sharp1 
            roundness[nstar] = around
            nstar = nstar+1

# endfor

# NSTAR is now the index of last star found
    nstar = nstar-1        

    print ' No. of sources rejected by SHARPNESS criteria'+badsharp
    print ' No. of sources rejected by ROUNDNESS criteria'+badround
    print ' No. of sources rejected by CENTROID  criteria'+badcntrd

# Any stars found?
    if nstar > 0:
        x = x[0:nstar]  
        y = y[0:nstar]
        flux = flux[0:nstar]
        sharp = sharp[0:nstar]  
        roundness = roundness[0:nstar]

        print ' Threshold above background for this pass was',hmin

        return([x, y, flux, sharp, roundness])                                      

    