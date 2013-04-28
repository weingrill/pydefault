'''
Created on 17.04.2013

@author: jwe
'''
class Spectrum(object):
    def __init__(self):
        import numpy as np
        mu, sigma = 1.0, 0.01 # mean and standard deviation
        r = 1000
        self.intensity = np.random.normal(mu, sigma, r)
        self.wavelength = np.linspace(350.0, 700.0, r)
        x = self.wavelength
        loc = 550.0
        scale = 3.0
        g = 0.2*np.exp( - (x - loc)**2 / (2 * scale**2) )
        self.intensity -= g
        loc = 656.0
        scale = 3.0
        g = 0.1*np.exp( - (x - loc)**2 / (2 * scale**2) )
        self.intensity += g
        self.cont = None
        
    def fromfile(self, filename):
        import pyfits
    
        filename = '/Users/jwe/Dropbox/Public/ses1/science20130408B-0030_botzfxsEcd.fits'
        hdulist = pyfits.open(filename)
        hdr = hdulist[0].header
        stars = hdulist[1].data
        hdulist.close()
    
    def continuum(self, width):   
        from numpy import array, ones, hstack
        from scipy.stats import hmean
        c = [] 
        leftpad = ones(width/2)*self.intensity[0]
        rightpad = ones(width-width/2)*self.intensity[-1]
        padded = hstack((leftpad, self.intensity, rightpad))

        for i in range(len(self.intensity)):
            p = padded[i:i+width]
            hm = hmean(p)
            c.append(hm*width)
        self.cont = array(c)
        
    def plot(self):
        import matplotlib

        matplotlib.use('WXAgg')
        import matplotlib.pyplot as plt
        plt.scatter(self.wavelength,self.intensity)
        plt.plot(self.wavelength,self.cont,'r')
        plt.show()
        
if __name__ == '__main__':
    

    f = '/Users/jwe/Dropbox/Public/ses1/science20130408B-0030_botzfxsEcd.fits'
    sp = Spectrum()
    sp.continuum(200)
    sp.plot()