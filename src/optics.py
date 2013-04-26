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

class detector(object):
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
            
    def width(self):
        """returns the width of the detectors size"""
        return self.pixels[0]*self.pixelsize

    def height(self):
        """returns the height of the detectors size"""
        return self.pixels[1]*self.pixelsize
    
    def size(self):
        """returns the dimensions of the sensor"""
        return [self.pixels[0]*self.pixelsize, self.pixels[1]*self.pixelsize]

    def diameter(self):
        """diagonal dimension of the sensor"""
        from math import sqrt
        return sqrt(self.pixels[0]**2+self.pixels[1]**2)*self.pixelsize 

    def fieldofview(self, telescope):
        """sensors field of view in radians for the diagonal, x and y"""
        from math import atan
        return (2.*atan(self.diameter()/(2.*telescope.focallength)),
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
    
class telescope(object):
    """
    telescope class
    """

    def __init__(self, diameter, focallength, name=None):
        """
        Constructor
        diameter and focal length in meters
        """
        teldb = {'C11':(0.2794,2.794), 
                 'MTO1000':(0.1,1.0)}
        if name in teldb:
            diameter, focallength = teldb[name]
        self.diameter = diameter
        self.focallength = focallength
    
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
        """retursn the angular field of view of a detector on the telescope"""
        from math import atan
        return 2.*atan(detector.diameter/(2.*self.focallength))
        
    def resolution(self, wavelength=550e-9):
        """returns the angular resolution in radians"""
        return 1.02*wavelength/self.diameter
    