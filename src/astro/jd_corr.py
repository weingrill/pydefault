#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Nov 21, 2016

@author: Leo Huckvale <leo.huckvale@postgrad.manchester.ac...> 

https://mail.scipy.org/pipermail/astropy/2014-April/003148.html

Hi all,

  Following recent discussion about implementing HJD on astropy-dev and 
github, I thought I might offer up one solution to correct timestamps to 
HJD/BJD. It seemed more appropriate here than on the dev mailing list.

  The code below defines the function "jd_corr" which takes Modified 
Julian Dates (on UTC scale), source RA and Dec and a correction type 
"bjd" or "hjd". This returns "new_jd", which is the corrected set of 
time stamps in an astropy.time.Time object. The delta-t corrections are 
obtained from barycentric ephemerides via jplephem, as recommended in 
the discussion on github, and a corresponding data package (in this 
case, de423).

Hopefully this might help someone. I'd appreciate any feedback - 
especially if I've overlooked something vital, which is entirely 
possible! I have relied in a large part on astropy.time to work its 
magic on time-standard conversions.

Best regards,
Leo
'''
import numpy as np

import astropy.time as time
import astropy.coordinates as coords
import astropy.units as u
import astropy.constants as const

import jplephem  # @UnresolvedImport
import de423  # @UnresolvedImport

VISTA_LAT = -24.6157    # Geographic coordinates of observatory
VISTA_LON = -70.3976

def jd_corr(mjd, ra, dec, jd_type='bjd'):
    """Return BJD or HJD for input MJD(UTC)."""
    
    # Initialise ephemeris from jplephem
    eph = jplephem.Ephemeris(de423)
    
    # Source unit-vector
    ## Assume coordinates in ICRS
    ## Set distance to unit (kilometers)
    src_vec = coords.ICRS(ra=ra, dec=dec,
                          unit=(u.degree, u.degree),
                          distance=coords.Distance(1, u.km))
    
    # Convert epochs to astropy.time.Time
    ## Assume MJD(UTC)
    t = time.Time(mjd, scale='utc', format='mjd', lat=VISTA_LAT, lon=VISTA_LON)  # @UndefinedVariable
    
    # Get Earth-Moon barycenter position
    ## NB: jplephem uses Barycentric Dynamical Time, e.g. JD(TDB)
    ## and gives positions relative to solar system barycenter
    barycenter_earthmoon = eph.position('earthmoon', t.tdb.jd)
    
    # Get Moon position vectors
    moonvector = eph.position('moon', t.tdb.jd)
    
    # Compute Earth position vectors
    pos_earth = (barycenter_earthmoon - moonvector * eph.earth_share)*u.km
    
    if jd_type == 'bjd':
        # Compute BJD correction
        ## Assume source vectors parallel at Earth and Solar System Barycenter
        ## i.e. source is at infinity
        corr = np.dot(pos_earth.T, src_vec.cartesian.value)/const.c
    elif jd_type == 'hjd':
        # Compute HJD correction via Sun ephemeris
        pos_sun = eph.position('sun', t.tdb.jd)*u.km
        sun_earth_vec = pos_earth - pos_sun
        corr = np.dot(sun_earth_vec.T, src_vec.cartesian.value)/const.c
    
    # TDB is the appropriate time scale for these ephemerides
    dt = time.TimeDelta(corr, scale='tdb', format='jd')  # @UndefinedVariable
    
    # Compute and return HJD/BJD as astropy.time.Time
    new_jd = t + dt
    
    return new_jd
