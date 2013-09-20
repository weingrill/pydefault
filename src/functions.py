'''
Created on Apr 26, 2013

@author: jwe
'''
def gaussian(x, mu = 0.0, sigma = 1.0):
    """ definition of the guassian function"""
    from math import sqrt, pi, exp
    y = 1 / ( sqrt(2. * pi * sigma) )  * exp( ((x - mu) / sigma)**2 )
    return y

def scaleto(values, bounds, k=None, d=None):
    """scales values within bounds"""
    import numpy as np
    x1, x2 = np.min(values),np.max(values)
    y1, y2 = bounds[0], bounds[1]
    if k is None:
        k =  (y2 - y1)/(x2 - x1)
    if d is None:
        d = y1 - k*x1
    return k*values + d

def histeq(im, nbr_bins=256):
    """http://www.janeriksolem.net/2009/06/histogram-equalization-with-python-and.html"""
    #get image histogram
    imhist,bins = histogram(im.flatten(),nbr_bins,normed=True)
    cdf = imhist.cumsum() #cumulative distribution function
    cdf = 255 * cdf / cdf[-1] #normalize

    #use linear interpolation of cdf to find new pixel values
    im2 = interp(im.flatten(),bins[:-1],cdf)

    return im2.reshape(im.shape), cdf