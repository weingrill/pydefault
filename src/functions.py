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

def normalize(values):
    """scales values to [0.0, 1.0]"""
    import numpy as np
    
    result = values - np.min(values)
    m = np.max(result)
    if abs(m) > 0.0:
        result /= m
    return result

def histeq(im, nbr_bins=256):
    """http://www.janeriksolem.net/2009/06/histogram-equalization-with-python-and.html"""
    from numpy import histogram, interp
    
    #get image histogram
    imhist,bins = histogram(im.flatten(),nbr_bins,normed=True)
    cdf = imhist.cumsum() #cumulative distribution function
    cdf = 255 * cdf / cdf[-1] #normalize

    #use linear interpolation of cdf to find new pixel values
    im2 = interp(im.flatten(),bins[:-1],cdf)

    return im2.reshape(im.shape), cdf

def sign(x):
    """returns the sign of a given number
    http://stackoverflow.com/questions/1986152/
    or use numpy.sign
    """
    #TODO: raise depreciation warning
    from math import copysign
    return copysign(1, x)
    
def rms(values):
    """returns the root mean square of values"""
    from numpy import mean, sqrt
    
    return sqrt(mean(values**2))

def modulus(image):
    """returns the modulus of an image,
    which can be used as an estimate for the sky background"""
    from numpy import histogram, argmax
    
    h, loc = histogram(image, bins=65536, range=[0,65535])
    m = argmax(h)
    return loc[m]

def logspace(lower, upper, num=20):
    """returns true logarithmic spacing between lower and upper
    instead of numpy version"""
    import numpy as np
    result = np.logspace(0.0, 1.0, num)
    return scaleto(result, [lower, upper])

def file_exists(filepath):
    ''' Check if a file exists and is accessible. '''
    try:
        f = open(filepath, 'rb')
    except IOError as e:
        return False
    f.close()
    return True