#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Joerg Weingrill"
__copyright__ = "Copyright 2018, Leibniz Institute for Astrophysics Potsdam (AIP)"
__credits__ = ["Joerg Weingrill"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Joerg Weingrill"
__email__ = "jweingrill@aip.de"
__status__ = "Development"
__date__="2018-04-12"

import matplotlib.pyplot as plt
from optics import radtoarcsec
import numpy as np
def resolution(diameter, wavelength=550e-9):
    """
    returns the angular resolution in radians
    https://en.wikipedia.org/wiki/Angular_resolution
    """
    return 1.220*wavelength/diameter

def dawes(diameter):
        """
        Dawes' limit
        this is an empirical result!
        D originally in inches
        result in arcsec
        https://en.wikipedia.org/wiki/Angular_resolution
        """
        return 116/(diameter*1000)

def pixelfov(pixelsize, focallength):
    """returns the angular resolution in radians"""
    return 2.*np.arctan(pixelsize/(2.*focallength))

def rayleigh(focalratio, wavelength=550e-9):
    """
    returns the Reyleigh criterion
    result in meter
    https://en.wikipedia.org/wiki/Angular_resolution
    """
    return 1.21966989*wavelength*focalratio
    
def airy(focalratio, wavelength=550e-9):
    """
    returns the linear diameter of the Airy disk
    result in meter
    https://en.wikipedia.org/wiki/Angular_resolution
    """
    return 2.0*1.21966989*wavelength*focalratio

def fwhm(focalratio, wavelength=550e-9):
    """
    returns the linear FWHM of a perfect stellar image
    result in meter
    """
    return 1.02*wavelength*focalratio

    
diameters = np.arange(0.05,0.51,0.01)
plt.title('Telescope resolution')
plt.plot(diameters*100, radtoarcsec(resolution(diameters)),label='resolution')
plt.plot(diameters*100, dawes(diameters),label='Dawes'' limit')
plt.xlabel('telescope diameter in cm')
plt.ylabel('resolution in arcsec')
plt.xticks(np.arange(5,51,5))
plt.ylim(0.0,2.5)
plt.xlim(5,50.0)
plt.grid()
plt.legend(loc='best')
plt.minorticks_on()
plt.savefig('/home/jwe/Pictures/Figures/telescope resolution.pdf')
plt.close()

focallengths = np.arange(0.5,5.0,0.01)

plt.title('Angular resolution')
pixelsizes = [2.4e-6, 3.75e-6, 4.65e-6, 5.6e-6, 9e-6]
colors = ['r','y', 'm', 'b', 'g']
for pixelsize,color in zip(pixelsizes, colors):
    plt.plot(focallengths, radtoarcsec(pixelfov(pixelsize, focallengths)),label=u'%.2f µm'%(pixelsize*1e6), color=color)
plt.xlabel('focal length in m')
plt.ylabel('angular resolution in arcsec')
plt.ylim(0.0, 2.5)
plt.xlim(0.5, 5.0)
plt.xticks(np.arange(0.5, 5.5, 0.5), ['%.1f'% fl for fl in np.arange(0.5, 5.5, 0.5)])
plt.grid()
plt.legend(loc='best')
plt.minorticks_on()
plt.savefig('/home/jwe/Pictures/Figures/angular resolution.pdf')
plt.close()


focalratios = np.arange(2.5,30.1,0.1)
plt.title('Spatial resolution')
plt.semilogx(focalratios, rayleigh(focalratios)*1e6, label='Reyleigh criterion', color='lightblue')
plt.semilogx(focalratios, airy(focalratios)*1e6, label='Airy disk diameter', color='lightgreen')
plt.semilogx(focalratios, fwhm(focalratios)*1e6, label='FWHM stellar image', color='pink')
for pixelsize,color in zip(pixelsizes, colors):
    plt.axhline(pixelsize*1e6, linestyle='--', color=color, linewidth=1.0)
    plt.text(30.0, pixelsize*1e6, u' %.2f µm' % (pixelsize*1e6), va='center')
plt.xlabel('focal ratio')
plt.ylabel(u'resolution in µm')
focalratios = np.array([2.5, 3.0, 3.3, 5.6, 6.3, 10.0, 15.0, 20.0, 30.0])
plt.xticks(focalratios, ['1/%.1f'%f for f in focalratios], rotation=270)
plt.ylim(0.0, 20.0)
plt.xlim(2.5,30.0)
plt.grid()
plt.legend(loc='best')
#plt.minorticks_on()
plt.savefig('/home/jwe/Pictures/Figures/spatial resolution.pdf')
plt.close()

def criticalfocal(diameter, pixelsize, wavelength=550e-9):
    """
    returns critical focal length (p. 22 in AIP)
    the focal length of the telescope must be smaller or equal to this.
    """
    return diameter*pixelsize/(0.51*wavelength)

plt.title('Critical focal length')
for pixelsize,color in zip(pixelsizes, colors):
    plt.plot(diameters*100,criticalfocal(diameters, pixelsize),label=u'%.2f µm'%(pixelsize*1e6), color=color)
plt.xlabel('telescope diameter in cm')
plt.ylabel('maximum focal length in m')
plt.grid()
plt.legend(loc='best')
plt.minorticks_on()
plt.savefig('/home/jwe/Pictures/Figures/critical focal length.pdf')

