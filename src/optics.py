'''
Created on Jan 24, 2013

@author: jwe
'''

def radtoarcsec(radians):
    """converts radians to arcseconds"""
    from math import pi
    return radians*3600.*180./pi

def radtoarcmin(radians):
    """converts radians to arcminutes"""
    from math import pi
    return radians*60.*180./pi

def radtodeg(radians):
    """converts radians to degrees"""
    from math import pi
    return radians*180./pi

class Detector(object):
    """
    detector class
    """
    def __init__(self, pixels, pixelsize, detname=None):
        """
        constructor of the detector class
        requires pixels with two entries: [width, height]
        pixelsize (assumed to be quadratic)
        """
        detdb = {'ICX205AL':([1280,960],4.65e-6),
                 'ICX204AL':([1024,768],4.65e-6),
                 'ICX098BL': ([640,480],5.6e-6),
                 'ICX618ALA':([640,480],5.6e-6),
                 }
        if detname in detdb:
            pixels, pixelsize = detdb[detname]
        self.pixels = pixels
        self.pixelsize = pixelsize
    
    @property        
    def width(self):
        """returns the width of the detectors size"""
        return self.pixels[0]*self.pixelsize

    @property        
    def height(self):
        """returns the height of the detectors size"""
        return self.pixels[1]*self.pixelsize
    
    @property        
    def size(self):
        """returns the dimensions of the sensor"""
        return [self.pixels[0]*self.pixelsize, self.pixels[1]*self.pixelsize]

    @property        
    def diameter(self):
        """diagonal dimension of the sensor"""
        from math import sqrt
        return sqrt(self.pixels[0]**2+self.pixels[1]**2)*self.pixelsize 

    def fieldofview(self, telescope):
        """sensors field of view in radians for the diagonal, x and y"""
        from math import atan
        return (2.*atan(self.diameter/(2.*telescope.focallength)),
                2.*atan(self.pixels[0]*self.pixelsize/(2.*telescope.focallength)),
                2.*atan(self.pixels[1]*self.pixelsize/(2.*telescope.focallength)))

    def pixelfov(self, telescope):
        """returns the angular resolution in radians"""
        from math import atan
        return 2.*atan(self.pixelsize/(2.*telescope.focallength))
    
    def criticalfocal(self, telescope, wavelength=550e-9):
        """
        returns critical focal length (p. 22 in AIP)
        the focal length of the telescope must be smaller or equal to this.
        """
        return telescope.diameter*self.pixelsize/(0.51*wavelength)
    
    def focusdepth(self, telescope, wavelength=550e-9):
        """returns the depth of focus (AIP, p.58)"""
        n = telescope.focalratio()
        return max([self.pixelsize*n, 4.0*n**2*wavelength])
    
class Telescope(object):
    """
    telescope class
    """

    def __init__(self, diameter, focallength, name=None):
        """
        diameter and focal length in meters 
        an optional given name overrides the diameter and the focallength
        """
        teldb = {'C11':          (0.2794, 2.794), 
                 'MTO1000':      (0.1, 1.0),
                 'Cassegrain50': (0.5, 7.5),
                 'Zeiss110':     (0.11, 0.11*7.5),
                 'BMK':          (0.30, 0.75),
                 'ETX':          (0.08, 0.4)}
        if name in teldb:
            diameter, focallength = teldb[name]
        self.diameter = diameter
        self.focallength = focallength
    
    @property
    def focalratio(self):
        """returns the focal ratio of the telescope"""
        return self.focallength/self.diameter
    
    def fwhm(self, wavelength=550e-9):
        """returns the linear FWHM of a perfect stellar image"""
        return 1.02*wavelength*self.focallength/self.diameter
    
    def airy(self, wavelength=550e-9):
        """returns the linear diameter of the Airy disk"""
        return 2.44*wavelength*self.focalratio()
    
    def fieldofview(self, detector):
        """returns the angular field of view of a detector on the telescope"""
        from math import atan
        return 2.*atan(detector.diameter/(2.*self.focallength))
        
    def resolution(self, wavelength=550e-9):
        """returns the angular resolution in radians"""
        return 1.02*wavelength/self.diameter
    
    def magnification(self, minpupil=0.001, maxpupil=0.0075):
        """returns the minimum and the maximum magnifications"""
        return (self.diameter/maxpupil, self.diameter/minpupil) 
    
    def limitingmagnitude(self):
        from math import log10
        return 5*log10(self.diameter) + 2.7   
 
class Eyepiece(object):
    def __init__(self, focallength, fov = 55):
        """
        focal length in meters
        fov in degrees
        """
        self.focallength = focallength
        self.fov = fov
    
    def magnification(self, telescope):
        return telescope.focallength/self.focallength
    
    def exitpupil(self, telescope):
        return telescope.diameter/self.magnification(telescope)
    
    def truefov(self, telescope):
        return self.fov/self.magnification(telescope)

if __name__ == '__main__':
    c11 = Telescope(0.0, 0.0, name='C11')
    minmag, maxmag = c11.magnification(maxpupil = 0.0075)
    print 'C11:'
    print '%i ... %i' % (minmag,maxmag)
    print 'ideal: %i mm' % (1000.0*c11.focallength/60.0)
    print 1000.0*c11.focallength/minmag, 1000.0*c11.focallength/maxmag, 'mm'
    
    print '\nCassegrain 50:'
    zeiss50 = Telescope(0.0, 0.0, name='Cassegrain50')
    minmag,maxmag = zeiss50.magnification()
    print '%i fach ... %i fach' % (minmag,maxmag)
    print 'ideal: %i mm' % (1000.0*zeiss50.focallength/60.0)
    print 'Okularbrennweite: %.1f mm ... %.1f mm' % (1000.0*zeiss50.focallength/minmag, 1000.0*zeiss50.focallength/maxmag)
    
    print '\nMeade 8":'
    meade = Telescope(0.2, 1.2)
    minmag, maxmag = meade.magnification(maxpupil = 0.007)
    print '%i ... %i' % (minmag,maxmag)
    print 'ideal: %i mm' % (1000.0*meade.focallength/60.0)
    print 1000.0*meade.focallength/minmag, 1000.0*meade.focallength/maxmag, 'mm'
    
    print '\nBMK:'
    bmk = Telescope(0.3, 0.75, name='BMK')
    print '1/f: %f' % bmk.focalratio
    ccd = Detector([10560,10560],9e-6)
    print 'FOV: %f diag.' % radtodeg(bmk.fieldofview(ccd))
    print 'FOV: %f x %f' % (radtodeg(ccd.fieldofview(bmk)[1]),radtodeg(ccd.fieldofview(bmk)[2]))
    print 'BMK resolution: %f"' % radtoarcsec(bmk.resolution())
    print 'BMK resolution: %f"' % radtoarcsec(ccd.pixelfov(bmk))

    print '\nETX:'
    etx = Telescope(0.07, 0.35, name='ETX')
    minmag, maxmag = etx.magnification(maxpupil = 0.007)
    print '%i ... %i' % (minmag,maxmag)
    print 'ideal: %i mm' % (1000.0*etx.focallength/60.0)
    print 'Okularbrennweite: %.1f mm ... %.1f mm' % (1000.0*etx.focallength/minmag, 1000.0*etx.focallength/maxmag)
    