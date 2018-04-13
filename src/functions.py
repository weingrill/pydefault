'''
Created on Apr 26, 2013

@author: jwe
'''
import numpy as np
import warnings

def gaussian(x, amp = 1.0, mu = 0.0, sigma = 1.0):
    """ definition of the statistical Gaussian function:
    amp: amplitude
    sigma**2: variance
    mu: average
    """
    # see scipy.signal.gaussian
    y = amp * np.exp(-0.5*((x - mu) / sigma)**2 ) / (sigma * np.sqrt(2. * np.pi * sigma))
    return y

def gauss(x, a, x0, sigma):
    """gauss function to be used by gauss_fit"""
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def gauss_fit(x, y, amp=None, mean=None, sigma=None ):
    '''
    fits a gaussian to (x, y)
    adopted from 
    http://stackoverflow.com/questions/19206332/gaussian-fit-for-python
    '''
    from scipy.optimize import curve_fit
    
    n = len(x)
    if amp is None:
        amp = max(y)
    if mean is None:
        mean = sum(x*y)/n
    if sigma is None:
        sigma = sum(y*(x-mean)**2)/n
    popt,_ = curve_fit(gauss, x, y, p0=[amp, mean, sigma])
    popt[2] = abs(popt[2])
    return popt  

def scaleto(values, bounds, k=None, d=None):
    """scales values within bounds"""
    x1, x2 = np.min(values),np.max(values)
    y1, y2 = bounds[0], bounds[1]
    if k is None:
        k =  (y2 - y1)/(x2 - x1)
    if d is None:
        d = y1 - k*x1
    return k*values + d

def normalize(values):
    """scales values to [0.0, 1.0]"""
    result = values - np.min(values)
    m = np.max(result)
    if abs(m) > 0.0:
        result /= m
    return result

def histeq(im, nbr_bins=256):
    """http://www.janeriksolem.net/2009/06/histogram-equalization-with-python-and.html"""
    
    #get image histogram
    imhist,bins = np.histogram(im.flatten(),nbr_bins,normed=True)
    cdf = imhist.cumsum() #cumulative distribution function
    cdf = 255 * cdf / cdf[-1] #normalize

    #use linear interpolation of cdf to find new pixel values
    im2 = np.interp(im.flatten(),bins[:-1],cdf)

    return im2.reshape(im.shape), cdf

def sign(x):
    """returns the sign of a given number
    http://stackoverflow.com/questions/1986152/
    or use numpy.sign
    """
    from math import copysign
    warnings.warn("Deprecated! please use numpy.sign")
    return copysign(1, x)
    
def rms(values):
    """returns the root mean square of values"""
    return np.sqrt(np.mean(values**2))

def modulus(image):
    """returns the modulus of an image,
    which can be used as an estimate for the sky background"""
    
    h, loc = np.histogram(image, bins=65536, range=[0,65535])
    m = np.argmax(h)
    return loc[m]

def logspace(lower, upper, num=20):
    """returns true logarithmic spacing between lower and upper
    instead of numpy version
    superseded by numpy.geomspace numpy>=v1.13"""
    warnings.warn("superseded by numpy.geomspace numpy>=v1.13")
    result = np.logspace(0.0, 1.0, num)
    return scaleto(result, [lower, upper])

def file_exists(filepath):
    ''' Check if a file exists and is accessible. '''
    warnings.warn("Deprecated! please use os.path.exists")
    try:
        f = open(filepath, 'rb')
    except IOError:
        return False
    f.close()
    return True

def autocorrelate(t, y):
    """
    calculate the autocorrelation on irregular data
    """
    from scipy.interpolate import interp1d

    #determine the median time step dt:
    dt = np.median(t-np.roll(t, 1))
    t0 = min(t)
    t1 = max(t)
    n = int((t1-t0)/dt)
    nt = np.linspace(t0, t1, n)
    f = interp1d(t, y, kind='linear')
    fn = f(nt)
    
    m = np.mean(fn)
    # we do padding
    k = 2
    x = fn - m
    # complex FFT, since we need the conjugate
    s = np.fft.fft(x, k*n)
    ac = np.real(np.fft.ifft(s*s.conjugate()))
    ac /= ac[0]
    lag = nt-nt[0]
    return ac[:n], lag[:n]

def sigma_clip(t, y, e=None,sigmas=3.0):
    """
    performs sigma clipping on a lightcurve
    """
    m = np.nanmean(y)
    s = np.nanstd(y)
    valid = abs(y-m)<sigmas*s
    t_clipped = np.compress(valid, t)
    y_clipped = np.compress(valid, y)
    if not e is None:
        e_clipped = np.compress(valid, e)
        return t_clipped, y_clipped, e_clipped
    else:
        return t_clipped, y_clipped

def phase(t, y, period):
    """
    returns the phased lightcurve
    """
    tp = t % period
    i = tp.argsort()
    return tp[i], y[i]

def smooth(x, n=101, width=2.0):
    """
    smooth using a gaussian kernel
    """
    from scipy import signal
    kernel = signal.gaussian(n, width)
    n = len(x)/2
    padded_x = np.pad(x, n, mode='mean')
    smoothed = signal.fftconvolve(padded_x, kernel, mode='same')[n:-n]/5
    return smoothed

def largmin(array, i0, direction='left'):
    """
    finds the local minimum based on direction
    """
    if direction not in ['left', 'right']:
        raise ValueError

    if i0 < 0 or i0 > len(array):
        raise ValueError('initial value out of bounds')
    i = i0
    if direction == 'left':
        while i > 0 and array[i-1] < array[i]:
            i -= 1
    elif direction == 'right':
        while i+1 < len(array) and array[i+1] < array[i]:
            i += 1
    return i

def runningmean(x, N=97):
    return np.convolve(x, np.ones((N,))/N, mode='same')
    
def regression(x, y, x_err, y_err):
    '''
    source 
    https://stackoverflow.com/questions/22670057/linear-fitting-in-python-with-uncertainty-in-both-x-and-y-coordinates
    '''
    from scipy.odr import Model, RealData, ODR
    # Initiate some data, giving some randomness using random.random().
    x = np.array(x)
    y = np.array(y)
    
    x_err = np.array(x_err)
    y_err = np.array(y_err)
    
    # Define a function (quadratic in our case) to fit the data with.
    def lin_func(p, x):
        k, d = p
        return k*x + d
    
    # Create a model for fitting.
    quad_model = Model(lin_func)
    
    # Create a RealData object using our initiated data from above.
    data = RealData(x, y, sx=x_err, sy=y_err)
    
    # Set up ODR with the model and data.
    odr = ODR(data, quad_model, beta0=[0., 1.])
    
    # Run the regression.
    out = odr.run()
    
    # Use the in-built pprint method to give us results.
    out.pprint()
    '''Beta: [ 1.01781493  0.48498006]
    Beta Std Error: [ 0.00390799  0.03660941]
    Beta Covariance: [[ 0.00241322 -0.01420883]
     [-0.01420883  0.21177597]]
    Residual Variance: 0.00632861634898189
    Inverse Condition #: 0.4195196193536024
    Reason(s) for Halting:
      Sum of squares convergence'''
    
    x_fit = np.linspace(x[0], x[-1], 1000)
    y_fit = lin_func(out.beta, x_fit)
    
    #plt.errorbar(x, y, xerr=x_err, yerr=y_err, linestyle='None', marker='x')
    #plt.plot(x_fit, y_fit)
    
    #plt.show()        
    return out.beta



