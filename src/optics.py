'''
Created on Jan 24, 2013

@author: jwe
'''

def radtoarcsec(radians):
    from math import pi
    return radians*3600.*180./pi

def radtoarcmin(radians):
    from math import pi
    return radians*60.*180./pi

def radtodeg(radians):
    from math import pi
    return radians*180./pi

class detector(object):
    """
    detector class
    """
    def __init__(self, pixels, pixelsize):
        """
        constructor of the detector class
        requires pixels with two entries: [width, height]
        pixelsize (assumed to be quadratic)
        """
        self.pixels = pixels
        self.pixelsize = pixelsize
    
    def width(self):
        return self.pixels[0]*self.pixelsize

    def height(self):
        return self.pixels[1]*self.pixelsize
    
    def size(self):
        return self.pixels*self.pixelsize

    def diameter(self):
        from math import sqrt
        return sqrt(sum(self.pixels**2))*self.pixelsize 

    def fieldofview(self, telescope):
        from math import atan
        return 2.*atan(self.diameter()/(2.*telescope.focallength))

    def pixelfov(self, telescope):
        """
        returns the angular 
        """
        from math import atan
        return 2.*atan(self.pixelsize/(2.*telescope.focallength))
    def criticalfocal(self, telescope, wavelength=550e-9):
        """
        returns critical focal length (p. 22 in AIP)
        the focal length of the telescope must be smaller or equal to this.
        """
        return telescope.diameter*self.pixelsize/(0.51*wavelength)
    
class telescope(object):
    """
    telescope class
    """

    def __init__(self, diameter, focallength):
        """
        Constructor
        diameter and focal length in meters
        """
        self.diameter = diameter
        self.focallength = focallength
    
    def focalratio(self):
        """
        returns the focal ratio of the telescope
        """
        return self.focallength/self.diameter
    
    def fwhm(self, wavelength=550e-9):
        """
        returns the linear FWHM of a perfect stellar image
        """
        return 1.02*wavelength*self.focallength/self.diameter
    
    def airy(self, wavelength=550e-9):
        """
        returns the linear diameter of the Airy disk
        """
        return 2.44*wavelength*self.focalratio()
    
    def fieldofview(self, detector):
        from math import atan
        return 2.*atan(detector.diameter/(2.*self.focallength))
        
    def resolution(self, wavelength=550e-9):
        """
        returns the angular resolution in radians
        """
        return 1.02*wavelength/self.diameter
    