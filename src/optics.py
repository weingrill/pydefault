#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jan 24, 2013

@author: joerg@weingrill.net
2013 - 2015 © Jörg Weingrill
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

class Angle(object):
    def __init__(self, value):
        """
        takes the value in radians
        """
        self.value = value
        self.precision = 2
    def __getitem__(self):
        return self.value
    
    def __str__(self):
        from math import pi
        result = self.value*180./pi
        if abs(result)>1.0:
            return u"%.2f˚" % result
        if abs(result*60.0)>1.0:
            return u"%.2fʹ" % (result*60.0)
        if abs(result*3600.0)>1.0:
            return u"%.2fʺ" % (result*3600.0)
        

class Detector(object):
    """
    detector class
    """
    def __init__(self, pixels=[], pixelsize=None, name=None, telescope=None):
        """
        constructor of the detector class
        requires pixels with two entries: [width, height]
        pixelsize (assumed to be quadratic)
        """
        detectors = {'ICX205AL':([1280,960],4.65e-6),
                 'ICX204AL':([1024,768],4.65e-6),
                 'ICX098BL': ([640,480],5.6e-6),
                 'ICX618ALA':([640,480],5.6e-6),
                 'AR0130CS':([1280,960],3.75e-6), #ASI120MM
                 'STA1600LN': ([10560,10560], 9e-6)
                 }
        if name in detectors:
            pixels, pixelsize = detectors[name]
        self.pixels = pixels
        self.pixelsize = pixelsize
        self.name = name
        self.telescope=telescope
    
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
        return [self.width, self.height]

    @property        
    def diagonal(self):
        """diagonal dimension of the sensor"""
        from math import sqrt
        return sqrt(self.pixels[0]**2+self.pixels[1]**2)*self.pixelsize 

    @property        
    def fieldofview(self):
        """sensors field of view in radians for the diagonal, x and y"""
        from math import atan
        return (2.*atan(self.diagonal/(2.*self.telescope.focallength)),
                2.*atan(self.pixels[0]*self.pixelsize/(2.*self.telescope.focallength)),
                2.*atan(self.pixels[1]*self.pixelsize/(2.*self.telescope.focallength)))

    @property        
    def pixelfov(self):
        """returns the angular resolution in radians"""
        from math import atan
        return 2.*atan(self.pixelsize/(2.*self.telescope.focallength))
    
    def criticalfocal(self, wavelength=550e-9):
        """
        returns critical focal length (p. 22 in AIP)
        the focal length of the telescope must be smaller or equal to this.
        """
        result = self.telescope.diameter*self.pixelsize/(0.51*wavelength) 
        assert(self.telescope.focallength <= result)
        return result
    
    def focusdepth(self, wavelength=550e-9):
        """returns the depth of focus (AIP, p.58)"""
        n = self.telescope.focalratio()
        return max([self.pixelsize*n, 4.0*n**2*wavelength])
    
    def __str__(self):
        data = {'name': self.name,
                'telname': self.telescope.name,
                'diagonal': self.diagonal * 1000,
                'fovx': radtodeg(self.fieldofview[1]),
                'fovy': radtodeg(self.fieldofview[2]),
                'pixelfov': radtoarcsec(self.pixelfov)
                }
        result = """name: %(name)s @ %(telname)s
diagonal:  %(diagonal).3f mm
FOV:       %(fovx).2f° × %(fovy).2f°
pixel FOV: %(pixelfov).3fʺ
        """ % data
        return result
    
class Telescope(object):
    """
    telescope class
    """

    def __init__(self, diameter=None, focallength=None, name=None):
        """
        diameter and focal length in meters 
        an optional given name overrides the diameter and the focallength
        """
        teldb = {'Celestron11':  (0.2794, 2.794), 
                 'MTO1000':      (0.1, 1.0),
                 'Cassegrain50': (0.5, 7.5),
                 'Zeiss110':     (0.11, 0.11*7.5),
                 'BMK':          (0.30, 0.75),
                 'ETX':          (0.08, 0.4)}
        if name in teldb:
            diameter, focallength = teldb[name]
        self.diameter = diameter
        self.focallength = focallength
        self.name = name
    
    @property
    def focalratio(self):
        """
        returns the focal ratio of the telescope
        """
        return self.focallength/self.diameter
    
    def fwhm(self, wavelength=550e-9):
        """
        returns the linear FWHM of a perfect stellar image
        result in meter
        """
        return 1.02*wavelength*self.focalratio
    
    @property
    def dawes(self):
        """
        Dawes' limit
        this is an empirical result!
        D originally in inches
        result in arcsec
        https://en.wikipedia.org/wiki/Angular_resolution
        """
        return 116/(self.diameter*1000)
    
    def rayleigh(self, wavelength=550e-9):
        """
        returns the Reyleigh criterion
        result in meter
        https://en.wikipedia.org/wiki/Angular_resolution
        """
        return 1.220*wavelength*self.focalratio
    
    def airy(self, wavelength=550e-9):
        """
        returns the linear diameter of the Airy disk
        result in meter
        https://en.wikipedia.org/wiki/Angular_resolution
        """
        return 2.44*wavelength*self.focalratio
    
    def fieldofview(self, detector):
        """
        returns the angular field of view of a detector on the telescope
        result in radians
        """
        from math import atan
        return 2.*atan(detector.diameter/(2.*self.focallength))
        
    def resolution(self, wavelength=550e-9):
        """
        returns the angular resolution in radians
        https://en.wikipedia.org/wiki/Angular_resolution
        """
        return 1.220*wavelength/self.diameter
    
    def magnification(self, minpupil=0.001, maxpupil=0.0075):
        """returns the minimum and the maximum magnifications"""
        return (self.diameter/maxpupil, self.diameter/minpupil) 
    
    @property
    def limitingmagnitude(self):
        from math import log10
        return 5*log10(self.diameter) + 2.7  
    
    def __str__(self): 
        minmag,maxmag = self.magnification()
        data = {'name': self.name,
                'focalratio': self.focalratio,
                'fwhm': self.fwhm()*1e6,
                'airy': self.airy()*1e6,
                'dawes': self.dawes,
                'resolution': radtoarcsec(self.resolution()),
                'minmag': minmag, 
                'maxmag': maxmag,
                'ideal': 1000.0*self.focallength/60.0, # 1arcsec seeing meets 1arcmin resolution of the eye
                'eyemax': 1000.0*self.focallength/minmag,
                'eyemin': 1000.0*self.focallength/maxmag}
        result =  """%(name)s:
focal ratio 1/f:    %(focalratio).1f
stellar fwhm:       %(fwhm).1f μm
Airy disc:          %(airy).1f μm
resolution:         %(resolution).3f ʺ
magnifications:     %(minmag)i × ... %(maxmag)i ×
ideal eyepiece:     %(ideal)i mm
range for eyepiece: %(eyemin).1f mm ... %(eyemax).1f mm""" % data
        
        return result

 
class Eyepiece(object):
    def __init__(self, focallength = None, 
                 fov = 55, 
                 name = None, 
                 telescope = None,
                 fieldstop = None):
        """
        Vergrößerung = Teleskopbrennweite : Okularbrennweite
        
        Austrittspupille (AP)
        AP = Teleskopöffnung : Vergrößerung
        Oder
        AP = Okularbrennweite : Öffnungszahl
        (Öffnungszahl eines f/10 Teleskopes ist 10)
        
        Wahres Gesichtsfeld:
        W = 2*Arctan(D/(2*F)
        D = Durchmesser der Feldblende
        F = Brennweite des Teleskops
        Oder (Faustformel)
        W = scheinbares Gesichtsfeld : Vergrößerung

        focal length in meters
        fov in degrees
        """
        eyepieces = {'TS ED 16mm': (0.016, 60, None),
                     'TS ED 15mm': (0.015, 60, None),
                     'TS ED 25mm': (0.025, 50, None),
                     'FMC 38mm':   (0.038, 70, None),
                     'Plössl 9mm': (0.009, 50, 0.007),
                     'ES WA 10mm': (0.010, 50, None)
                     }
        if name in eyepieces:
            focallength, fov, fieldstop = eyepieces[name]
        else:
            print eyepieces.keys()
            raise ValueError
        self.focallength = focallength
        self.fov = fov
        self.telescope = telescope
        self.fieldstop = fieldstop
        self.name = name
    
    @property
    def magnification(self):
        """
        returns the magnification in times (×)
        """
        return self.telescope.focallength/self.focallength
    
    @property
    def exitpupil(self):
        """
        returns the exitpupil in meters
        """
        return self.telescope.diameter/self.magnification
    
    @property
    def truefov(self):
        """
        returns the true field of view in degrees
        """
        from math import atan
        if self.fieldstop is None:
            return self.fov/self.magnification
        else:
            return radtodeg(2*atan(self.fieldstop/(2*self.telescope.focallength)))
    
    def __str__(self):
        data = {'name': self.name,
                'telname': self.telescope.name,
                'exitpupil': self.exitpupil*1000.0,
                'magnification': self.magnification,
                'truefov': self.truefov*60.0}
        result = """%(name)s @ %(telname)s:
exitpupil:     %(exitpupil).1f mm
magnification: %(magnification)d× 
true FOV:      %(truefov)d '
        """ % data
        return result
        
    

    