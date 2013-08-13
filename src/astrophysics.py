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

r_sun = 6.96342e8
m_sun = 1.9891e30

class Star(object):
    def __init__(self, mass = 1.9891e30):
        """creates a star with given mass in kg. If mass is less then 50 then
        solar masses as units are assumed"""
        self.mass = mass # kg
        if mass<=50:
            mass *=m_sun
            self.mass = mass
        #if mass<=2*m_sun:
        #    self.radius = r_sun*(mass/m_sun)**0.9
        #if mass>2*m_sun and mass<20*m_sun:
        #    self.radius = r_sun*(mass/m_sun)**(15/19)
    def get_radius(self):
        if self.mass<=2*m_sun:
            self._radius = r_sun*(self.mass/m_sun)**0.9
        if self.mass>2*m_sun and self.mass<20*m_sun:
            self._radius = r_sun*(self.mass/m_sun)**(15/19)
        return self._radius
    
    #@radius.setter
    def set_radius(self, value):
        self._radius = value
        self.mass = m_sun*(self._radius/r_sun)**(19/15)
        
    radius = property(get_radius, set_radius)
        
    def luminosity(self):
        """returns the luminosity in Watts
        taken from:
        http://physics.ucsd.edu/students/courses/winter2008/managed/physics223/documents/Lecture7%13Part3.pdf
        """
        l_sun = 3.846e26 # [W] 
        if self.mass>=2*m_sun and self.mass<20*m_sun:
            return l_sun*(self.mass/m_sun)**3
        
        if self.mass>=0.5*m_sun and self.mass<2*m_sun:
            return l_sun*(self.mass/m_sun)**5.5/(self.radius/r_sun)**0.5
    
    def surface_gravity(self):
        """surface gravity of the star in meters per socond squared"""
        f = 6.673e-11 # gravitational constant [m**3 kg**-1 s**-2]
        return f*self.mass/self.radius**2 # [m s**-2]
    
    def breakup(self):
        """returns the breakup period in seconds"""
        from numpy import sqrt, pi
        f = 6.673e-11 # gravitational constant [m**3 kg**-1 s**-2]
        return sqrt(4 * pi**2 * self.radius**3/(f*self.mass))
    
    def lifetime(self):
        """returns the lifetime on the main-sequence in years"""
        m_sun = 1.9891e30
        return 10**10 * (self.mass/m_sun)**-2.5 #[years]
        
    def t_eff(self):
        """returns the effective Temperature in Kelvin"""
        from numpy import pi
        s = 5.6704e-8 #[W m**-2 K**4]
        return (self.luminosity()/(4*pi*self.radius**2*s))**0.25
    
    def spectral_class(self):
        """returns the spectral class of the star"""
        t = self.t_eff()
        sp_types = ['O','B','A','F','G','K','M','L','T','Y']
        t_lower = [33000,10000,7500,6000,5200,3700,2000,1300,700,0]
        spt = None
        for tl in t_lower:
            if t>=tl:
                i = t_lower.index(tl)
                spt = sp_types[i]
                break
        # subtype
        sub=0
        if spt in ['B','A','F','G','M','L','T']:
            sub = 10*(1.0-(t-t_lower[i])/(t_lower[i-1]-t_lower[i]))
        if spt=='K':
            sub = 7*(1.0-(t-t_lower[i])/(t_lower[i-1]-t_lower[i]))
            
        return '%s%i' % (spt,sub)
    
    def absolute_magnitude(self):
        from numpy import log10
        l_sun = 3.846e26
        return 4.83 -2.512*log10(self.luminosity()/l_sun)
    
if __name__ == '__main__':
    star = Star(1.05)
    print star.radius
    print '%.2f days' % (star.breakup()/86400.)
    print star.surface_gravity()
    print '%i K' % star.t_eff()
    print star.spectral_class()
    print star.lifetime()
    print star.absolute_magnitude()
    
    