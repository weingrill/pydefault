'''
Created on 17.04.2013

@author: jwe
'''
import numpy as np
class Spectrum(object):
    def __init__(self):
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
    
        filename = '/work1/jwe/Dropbox/Public/ses1/science20130410B-0008_botzfxsEcd.fits'
        hdulist = pyfits.open(filename)
        hdr = hdulist[0].header
        data = hdulist[0].data
        self.spectrum = data[0,40]
        self.raw = data[1,40]
        self.sigma = data[2,40]
        self.intensity = data[1,40]
        self.wavelength = np.arange(len(self.intensity))
        print self.wavelength.shape
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
        import matplotlib.pyplot as plt
        #plt.scatter(self.wavelength,self.intensity)
        plt.plot(self.spectrum,'r')
        plt.plot(self.raw,'b')
        plt.show()
        
if __name__ == '__main__':
    import pyspeckit
    from pylab import *
    #import wav2rgb
    
    speclist = pyspeckit.wrappers.load_IRAF_multispec('/work1/jwe/Dropbox/Public/ses1/science20130410B-0008_botzfxsEcd.fits')
    
    for spec in speclist:
        spec.units="Counts"
    
    SP = pyspeckit.Spectra(speclist)
    SPa = pyspeckit.Spectra(speclist,xunits='angstroms',quiet=False)
    
    SP.plotter(figure=figure(1))
    SPa.plotter(figure=figure(2))
    
    figure(3)
    clf()
    figure(4)
    clf()
    
    clr = [list(clr) for clr in matplotlib.cm.brg(linspace(0,1,len(speclist)))]
    #clr = [wav2rgb.wav2RGB(c) + [1.0] for c in linspace(380,780,len(speclist))][::-1]
    for ii,(color,spec) in enumerate(zip(clr,speclist)):
        spec.plotter(figure=figure(3), clear=False, reset=False, color=color, refresh=False)
    
        fig4=figure(4)
        fig4.subplots_adjust(hspace=0.35,top=0.97,bottom=0.03)
        spec.plotter(axis=subplot(10,1,ii%10+1), clear=False, reset=False, color=color, refresh=False)
        spec.plotter.axis.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(4) )
    
        if ii % 10 == 9:
            spec.plotter.refresh()
            spec.plotter.savefig('vega_subplots_%03i.png' % (ii/10+1))
            clf()
    
    spec.plotter.refresh()    

    f = '/Users/jwe/Dropbox/Public/ses1/science20130408B-0030_botzfxsEcd.fits'
    sp = Spectrum()
    #sp.continuum(200)
    sp.fromfile(f)
    sp.plot()