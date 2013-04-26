'''
Created on Apr 26, 2013

@author: jwe

'''

def color(Teff):
    """
    Cameron Reed of Alma College (Michigan) in 
    "The Composite Observational-Theoretical HR Diagram"
    The Journal of the Royal Astronomical Society of Canada
    February/March 1998 Volume92 Number 1 [669] page36
    has give an empirical fit of Mbol and CI to their Teff dependence.
    http://www.bogan.ca/astro/stars/star_lum.html
    """
    from math import log10
    logT = log10(Teff)
    if logT < 3.961:
        bv = -3.684*logT + 14.551
    else:
        bv = 0.344*logT**2 - 3.402*logT + 8.037
    return bv

def bc(Teff):
    """
    calculate the bolometric correction as a function of 
    effective temperature
    see 
      http://www.bogan.ca/astro/stars/star_lum.html
    for details
    """
    from math import log10
    logT4 = log10(Teff)-4.0
    return -8.499*logT4**4 + 13.421*logT4**3 - 8.131*logT4**2 - 3.901*logT4 - 0.438
    