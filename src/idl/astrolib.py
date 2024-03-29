#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Mar 30, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
def airtovac(wave_air):
    from numpy import where,array
    if type(wave_air) is float:
        wave_vac = wave_air
        for _ in range(2):
            sigma2 = (1e4/wave_vac)**2.
            fact = 1. +  5.792105e-2/(238.0185e0 - sigma2) + \
                                1.67917e-3/( 57.362e0 - sigma2)
            wave_vac = wave_air*fact              #Convert Wavelength
    elif type(wave_air) is list:
        
        wave_air = array(wave_air)
        wave_vac = array(wave_air)
    
        gi = where(wave_vac >= 2000)[0]     #Only modify above 2000 A
        for g in gi:
            sigma2 = (1e4/float(wave_vac[g]) )**2.     #Convert to wavenumber squared
            for _ in range(2):
                # Compute conversion factor
                fact = 1. +  5.792105e-2/(238.0185e0 - sigma2) + \
                                    1.67917e-3/( 57.362e0 - sigma2)
            
                wave_vac[g] = wave_air[g]*fact              #Convert Wavelength
    return wave_vac          

def premat(equinox1, equinox2, FK4 = True):
    """
    ;+
    ; NAME:
    ;       PREMAT
    ; PURPOSE:
    ;       Return the precession matrix needed to go from EQUINOX1 to EQUINOX2.  
    ; EXPLANTION:
    ;       This matrix is used by the procedures PRECESS and BARYVEL to precess 
    ;       astronomical coordinates
    ;
    ; CALLING SEQUENCE:
    ;       matrix = PREMAT( equinox1, equinox2, [ /FK4 ] )
    ;
    ; INPUTS:
    ;       EQUINOX1 - Original equinox of coordinates, numeric scalar.  
    ;       EQUINOX2 - Equinox of precessed coordinates.
    ;
    ; OUTPUT:
    ;      matrix - double precision 3 x 3 precession matrix, used to precess
    ;               equatorial rectangular coordinates
    ;
    ; OPTIONAL INPUT KEYWORDS:
    ;       /FK4   - If this keyword is set, the FK4 (B1950.0) system precession
    ;               angles are used to compute the precession matrix.   The 
    ;               default is to use FK5 (J2000.0) precession angles
    ;
    ; EXAMPLES:
    ;       Return the precession matrix from 1950.0 to 1975.0 in the FK4 system
    ;
    ;       IDL> matrix = PREMAT( 1950.0, 1975.0, /FK4)
    ;
    ; PROCEDURE:
    ;       FK4 constants from "Computational Spherical Astronomy" by Taff (1983), 
    ;       p. 24. (FK4). FK5 constants from "Astronomical Almanac Explanatory
    ;       Supplement 1992, page 104 Table 3.211.1.
    ;
    ; REVISION HISTORY
    ;       Written, Wayne Landsman, HSTX Corporation, June 1994
    ;       Converted to IDL V5.0   W. Landsman   September 1997
    ;-    
    """
    from numpy import pi,sin, cos, zeros
    deg_to_rad = pi/180.0
    sec_to_rad = deg_to_rad/3600.0
    
    t = 0.001e0*( equinox2 - equinox1)
    
    if not FK4:  
        st = 0.001e0*( equinox1 - 2000.)
        #  Compute 3 rotation angles
        a = sec_to_rad * t * (23062.181e0 + st*(139.656e0 +0.0139e0*st) \
            + t*(30.188e0 - 0.344e0*st+17.998e0*t))
        
        b = sec_to_rad * t * t * (79.280e0 + 0.410e0*st + 0.205e0*t) + a
        
        c = sec_to_rad * t * (20043.109e0 - st*(85.33e0 + 0.217e0*st) \
              + t*(-42.665e0 - 0.217e0*st -41.833e0*t))

    else:  
        st = 0.001e0*( equinox1 - 1900.e0)
        #  Compute 3 rotation angles
        
        a = sec_to_rad * t * (23042.53e0 + st*(139.75e0 +0.06e0*st) \
         + t*(30.23e0 - 0.27e0*st+18.0e0*t))
        
        b = sec_to_rad * t * t * (79.27e0 + 0.66e0*st + 0.32e0*t) + a
        
        c = sec_to_rad * t * (20046.85e0 - st*(85.33e0 + 0.37e0*st) \
           + t*(-42.67e0 - 0.37e0*st -41.8e0*t))


    sina, sinb,sinc = sin(a), sin(b), sin(c)
    cosa, cosb, cosc = cos(a), cos(b), cos(c)
    
    r = zeros(3,3)
    r[0,0] = [ cosa*cosb*cosc-sina*sinb, sina*cosb+cosa*sinb*cosc,  cosa*sinc]
    r[0,1] = [-cosa*sinb-sina*cosb*cosc, cosa*cosb-sina*sinb*cosc, -sina*sinc]
    r[0,2] = [-cosb*sinc, -sinb*sinc, cosc]
    
    return r

def baryvel( dje, deq, dvelh, dvelb, JPL = False):
    """
;+
; NAME:
;       BARYVEL
; PURPOSE:
;       Calculates heliocentric and barycentric velocity components of Earth.
;
; EXPLANATION:
;       BARYVEL takes into account the Earth-Moon motion, and is useful for 
;       radial velocity work to an accuracy of  ~1 m/s.
;
; CALLING SEQUENCE:
;       BARYVEL, dje, deq, dvelh, dvelb, [ JPL =  ] 
;
; INPUTS:
;       DJE - (scalar) Julian ephemeris date.
;       DEQ - (scalar) epoch of mean equinox of dvelh and dvelb. If deq=0
;               then deq is assumed to be equal to dje.
; OUTPUTS: 
;       DVELH: (vector(3)) heliocentric velocity component. in km/s 
;       DVELB: (vector(3)) barycentric velocity component. in km/s
;
;       The 3-vectors DVELH and DVELB are given in a right-handed coordinate 
;       system with the +X axis toward the Vernal Equinox, and +Z axis 
;       toward the celestial pole.      
;
; OPTIONAL KEYWORD SET:
;       JPL - if /JPL set, then BARYVEL will call the procedure JPLEPHINTERP
;             to compute the Earth velocity using the full JPL ephemeris.   
;             The JPL ephemeris FITS file JPLEPH.405 must exist in either the 
;             current directory, or in the directory specified by the 
;             environment variable ASTRO_DATA.   Alternatively, the JPL keyword
;             can be set to the full path and name of the ephemeris file.
;             A copy of the JPL ephemeris FITS file is available in
;                 http://idlastro.gsfc.nasa.gov/ftp/data/         
; PROCEDURES CALLED:
;       Function PREMAT() -- computes precession matrix
;       JPLEPHREAD, JPLEPHINTERP, TDB2TDT - if /JPL keyword is set
; NOTES:
;       Algorithm taken from FORTRAN program of Stumpff (1980, A&A Suppl, 41,1)
;       Stumpf claimed an accuracy of 42 cm/s for the velocity.    A 
;       comparison with the JPL FORTRAN planetary ephemeris program PLEPH
;       found agreement to within about 65 cm/s between 1986 and 1994
;
;       If /JPL is set (using JPLEPH.405 ephemeris file) then velocities are 
;       given in the ICRS system; otherwise in the FK4 system.   
; EXAMPLE:
;       Compute the radial velocity of the Earth toward Altair on 15-Feb-1994
;          using both the original Stumpf algorithm and the JPL ephemeris
;
;       IDL> jdcnv, 1994, 2, 15, 0, jd          ;==> JD = 2449398.5
;       IDL> baryvel, jd, 2000, vh, vb          ;Original algorithm
;               ==> vh = [-17.07243, -22.81121, -9.889315]  ;Heliocentric km/s
;               ==> vb = [-17.08083, -22.80471, -9.886582]  ;Barycentric km/s
;       IDL> baryvel, jd, 2000, vh, vb, /jpl   ;JPL ephemeris
;               ==> vh = [-17.07236, -22.81126, -9.889419]  ;Heliocentric km/s
;               ==> vb = [-17.08083, -22.80484, -9.886409]  ;Barycentric km/s
;
;       IDL> ra = ten(19,50,46.77)*15/!RADEG    ;RA  in radians
;       IDL> dec = ten(08,52,3.5)/!RADEG        ;Dec in radians
;       IDL> v = vb[0]*cos(dec)*cos(ra) + \   ;Project velocity toward star
;               vb[1]*cos(dec)*sin(ra) + vb[2]*sin(dec) 
;
; REVISION HISTORY:
;       Jeff Valenti,  U.C. Berkeley    Translated BARVEL.FOR to IDL.
;       W. Landsman, Cleaned up program sent by Chris McCarthy (SfSU) June 1994
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Added /JPL keyword  W. Landsman   July 2001
;       Documentation update W. Landsman Dec 2005
;-
    """
    from numpy import pi, reshape, cos, sin, sqrt
    # if N_params() LT 4 then begin
    #        print,'Syntax: BARYVEL, dje, deq, dvelh, dvelb'
    #        print,'    dje - input Julian ephemeris date'
    #        print,'    deq - input epoch of mean equinox of dvelh and dvelb'
    #        print,'    dvelh - output vector(3) heliocentric velocity comp in km/s' 
    #        print,'    dvelb - output vector(3) barycentric velocity comp in km/s'
    #        return
    # endif
    #Define constants
    dc2pi = 2*pi 
    #cc2pi = 2*pi 
    dc1 = 1.0E0
    dcto = 2415020.0e0
    dcjul = 36525.0e0                     # days in Julian year
    dcbes = 0.313e0
    dctrop = 365.24219572e0                # days in tropical year (...572 insig)
    dc1900 = 1900.0e0
    AU = 1.495978e8

    #Constants dcfel(i,k) of fast changing elements.
    dcfel = [1.7400353e00, 6.2833195099091e02,  5.2796e-6 \
          ,6.2565836e00, 6.2830194572674e02, -2.6180e-6 \
          ,4.7199666e00, 8.3997091449254e03, -1.9780e-5 \
          ,1.9636505e-1, 8.4334662911720e03, -5.6044e-5 \
          ,4.1547339e00, 5.2993466764997e01,  5.8845e-6 \
          ,4.6524223e00, 2.1354275911213e01,  5.6797e-6 \
          ,4.2620486e00, 7.5025342197656e00,  5.5317e-6 \
          ,1.4740694e00, 3.8377331909193e00,  5.6093e-6 ]
    dcfel = reshape(dcfel,[3,8])

    #constants dceps and ccsel(i,k) of slowly changing elements.
    dceps = [4.093198e-1, -2.271110e-4, -2.860401e-8 ]
    ccsel = [1.675104E-2, -4.179579E-5, -1.260516E-7 \
          ,2.220221E-1,  2.809917E-2,  1.852532E-5 \
          ,1.589963E00,  3.418075E-2,  1.430200E-5 \
          ,2.994089E00,  2.590824E-2,  4.155840E-6 \
          ,8.155457E-1,  2.486352E-2,  6.836840E-6 \
          ,1.735614E00,  1.763719E-2,  6.370440E-6 \
          ,1.968564E00,  1.524020E-2, -2.517152E-6 \
          ,1.282417E00,  8.703393E-3,  2.289292E-5 \
          ,2.280820E00,  1.918010E-2,  4.484520E-6 \
          ,4.833473E-2,  1.641773E-4, -4.654200E-7 \
          ,5.589232E-2, -3.455092E-4, -7.388560E-7 \
          ,4.634443E-2, -2.658234E-5,  7.757000E-8 \
          ,8.997041E-3,  6.329728E-6, -1.939256E-9 \
          ,2.284178E-2, -9.941590E-5,  6.787400E-8 \
          ,4.350267E-2, -6.839749E-5, -2.714956E-7 \
          ,1.348204E-2,  1.091504E-5,  6.903760E-7 \
          ,3.106570E-2, -1.665665E-4, -1.590188E-7 ]
    ccsel = reshape(ccsel,[3,17])

    #Constants of the arguments of the short-period perturbations.
    dcargs = [5.0974222e0, -7.8604195454652e2 \
           ,3.9584962e0, -5.7533848094674e2 \
           ,1.6338070e0, -1.1506769618935e3 \
           ,2.5487111e0, -3.9302097727326e2 \
           ,4.9255514e0, -5.8849265665348e2 \
           ,1.3363463e0, -5.5076098609303e2 \
           ,1.6072053e0, -5.2237501616674e2 \
           ,1.3629480e0, -1.1790629318198e3 \
           ,5.5657014e0, -1.0977134971135e3 \
           ,5.0708205e0, -1.5774000881978e2 \
           ,3.9318944e0,  5.2963464780000e1 \
           ,4.8989497e0,  3.9809289073258e1 \
           ,1.3097446e0,  7.7540959633708e1 \
           ,3.5147141e0,  7.9618578146517e1 \
           ,3.5413158e0, -5.4868336758022e2 ]
    dcargs = reshape(dcargs,[2,15])

    #Amplitudes ccamps(n,k) of the short-period perturbations.
    ccamps = \
    [-2.279594E-5,  1.407414E-5,  8.273188E-6,  1.340565E-5, -2.490817E-7 \
    ,-3.494537E-5,  2.860401E-7,  1.289448E-7,  1.627237E-5, -1.823138E-7 \
    , 6.593466E-7,  1.322572E-5,  9.258695E-6, -4.674248E-7, -3.646275E-7 \
    , 1.140767E-5, -2.049792E-5, -4.747930E-6, -2.638763E-6, -1.245408E-7 \
    , 9.516893E-6, -2.748894E-6, -1.319381E-6, -4.549908E-6, -1.864821E-7 \
    , 7.310990E-6, -1.924710E-6, -8.772849E-7, -3.334143E-6, -1.745256E-7 \
    ,-2.603449E-6,  7.359472E-6,  3.168357E-6,  1.119056E-6, -1.655307E-7 \
    ,-3.228859E-6,  1.308997E-7,  1.013137E-7,  2.403899E-6, -3.736225E-7 \
    , 3.442177E-7,  2.671323E-6,  1.832858E-6, -2.394688E-7, -3.478444E-7 \
    , 8.702406E-6, -8.421214E-6, -1.372341E-6, -1.455234E-6, -4.998479E-8 \
    ,-1.488378E-6, -1.251789E-5,  5.226868E-7, -2.049301E-7,  0.E0 \
    ,-8.043059E-6, -2.991300E-6,  1.473654E-7, -3.154542E-7,  0.E0 \
    , 3.699128E-6, -3.316126E-6,  2.901257E-7,  3.407826E-7,  0.E0 \
    , 2.550120E-6, -1.241123E-6,  9.901116E-8,  2.210482E-7,  0.E0 \
    ,-6.351059E-7,  2.341650E-6,  1.061492E-6,  2.878231E-7,  0.E0 ]
    ccamps = reshape(ccamps,[5,15])

    #Constants csec3 and ccsec(n,k) of the secular perturbations in longitude.
    ccsec3 = -7.757020E-8
    ccsec = [1.289600E-6, 5.550147E-1, 2.076942E00 \
          ,3.102810E-5, 4.035027E00, 3.525565E-1 \
          ,9.124190E-6, 9.990265E-1, 2.622706E00 \
          ,9.793240E-7, 5.508259E00, 1.559103E01 ]
    ccsec = reshape(ccsec,3,4)

    #Sidereal rates.
    dcsld = 1.990987e-7                   #sidereal rate in longitude
    ccsgd = 1.990969e-7                   #sidereal rate in mean anomaly

    #Constants used in the calculation of the lunar contribution.
    cckm = 3.122140e-5
    ccmld = 2.661699e-6
    ccfdi = 2.399485e-7

    #Constants dcargm(i,k) of the arguments of the perturbations of the motion
    # of the moon.
    dcargm = [5.1679830e0,  8.3286911095275e3 \
           ,5.4913150e0, -7.2140632838100e3 \
           ,5.9598530e0,  1.5542754389685e4 ]
    dcargm = reshape(dcargm,[2,3])

    #Amplitudes ccampm(n,k) of the perturbations of the moon.
    ccampm = [ 1.097594E-1, 2.896773E-7, 5.450474E-2,  1.438491E-7 \
           ,-2.223581E-2, 5.083103E-8, 1.002548E-2, -2.291823E-8 \
           , 1.148966E-2, 5.658888E-8, 8.249439E-3,  4.063015E-8 ]
    ccampm = reshape(ccampm, [4,3])

    #ccpamv(k)=a*m*dl,dt (planets), dc1mme=1-mass(earth+moon)
    ccpamv = [8.326827E-11, 1.843484E-11, 1.988712E-12, 1.881276E-12]
    dc1mme = 0.99999696e0

    #Time arguments.
    dt = (dje - dcto) / dcjul
    tvec = [1e0, dt, dt*dt]
    
    #Values of all elements for the instant(aneous?) dje.
    temp = (tvec.T * dcfel) % dc2pi
    dml = temp[0]
    forbel = temp[1:7]
    g = forbel[0]                         #old fortran equivalence
    
    deps = sum(tvec*dceps) % dc2pi
    sorbel = (tvec.T * ccsel) % dc2pi
    e = sorbel[0]                         #old fortran equivalence
    
    #Secular perturbations in longitude.
    dummy = cos(2.0)
    sn = sin((tvec[0:1].T * ccsec[1:2,:]) % dc2pi)

    #Periodic perturbations of the emb (earth-moon barycenter).
    pertl = sum(ccsec[0,:] * sn) + dt*ccsec3*sn[2]
    pertld = 0.0
    pertr = 0.0
    pertrd = 0.0
    for k in range(15):
        a = (dcargs[0,k]+dt*dcargs[1,k]) % dc2pi
        cosa = cos(a)
        sina = sin(a)
        pertl = pertl + ccamps[0,k]*cosa + ccamps[1,k]*sina
        pertr = pertr + ccamps[2,k]*cosa + ccamps[3,k]*sina
        if k < 11:
            pertld = pertld + (ccamps[1,k]*cosa-ccamps[0,k]*sina)*ccamps[4,k]
            pertrd = pertrd + (ccamps[3,k]*cosa-ccamps[2,k]*sina)*ccamps[4,k]

#Elliptic part of the motion of the emb.
    phi = (e*e/4e0)*(((8e0/e)-e)*sin(g) +5*sin(2*g) +(13/3e0)*e*sin(3*g))
    f = g + phi
    sinf = sin(f)
    cosf = cos(f)
    dpsi = (dc1 - e*e) / (dc1 + e*cosf)
    phid = 2*e*ccsgd*((1 + 1.5*e*e)*cosf + e*(1.25 - 0.5*sinf*sinf))
    psid = ccsgd*e*sinf / sqrt(dc1 - e*e)

#Perturbed heliocentric motion of the emb.
    d1pdro = dc1+pertr
    drd = d1pdro * (psid + dpsi*pertrd)
    drld = d1pdro*dpsi * (dcsld+phid+pertld)
    dtl = (dml + phi + pertl) % dc2pi
    dsinls = sin(dtl)
    dcosls = cos(dtl)
    dxhd = drd*dcosls - drld*dsinls
    dyhd = drd*dsinls + drld*dcosls

    #Influence of eccentricity, evection and variation on the geocentric
    # motion of the moon.
    pertl = 0.0
    pertld = 0.0
    pertp = 0.0
    pertpd = 0.0
    for k in [0,1,2]:
        a = (dcargm[0,k] + dt*dcargm[1,k]) % dc2pi
        sina = sin(a)
        cosa = cos(a)
        pertl = pertl + ccampm[0,k]*sina
        pertld = pertld + ccampm[1,k]*cosa
        pertp = pertp + ccampm[2,k]*cosa
        pertpd = pertpd - ccampm[3,k]*sina

    #Heliocentric motion of the earth.
    tl = forbel[1] + pertl
    sinlm = sin(tl)
    coslm = cos(tl)
    sigma = cckm / (1.0 + pertp)
    a = sigma*(ccmld + pertld)
    b = sigma*pertpd
    dxhd = dxhd + a*sinlm + b*coslm
    dyhd = dyhd - a*coslm + b*sinlm
    dzhd= -sigma*ccfdi*cos(forbel[2])

    #Barycentric motion of the earth.
    dxbd = dxhd*dc1mme
    dybd = dyhd*dc1mme
    dzbd = dzhd*dc1mme
    for k in range(3):
        plon = forbel[k+3]
        pomg = sorbel[k+1]
        pecc = sorbel[k+9]
        tl = (plon + 2.0*pecc*sin(plon-pomg)) % dc2pi
        dxbd = dxbd + ccpamv[k]*(sin(tl) + pecc*sin(pomg))
        dybd = dybd - ccpamv[k]*(cos(tl) + pecc*cos(pomg))
        dzbd = dzbd - ccpamv[k]*sorbel[k+13]*cos(plon - sorbel[k+5])

    #Transition to mean equator of date.
    dcosep = cos(deps)
    dsinep = sin(deps)
    dyahd = dcosep*dyhd - dsinep*dzhd
    dzahd = dsinep*dyhd + dcosep*dzhd
    dyabd = dcosep*dybd - dsinep*dzbd
    dzabd = dsinep*dybd + dcosep*dzbd

    #Epoch of mean equinox (deq) of zero implies that we should use
    # Julian ephemeris date (dje) as epoch of mean equinox.
    if deq == 0:
        dvelh = AU * [dxhd, dyahd, dzahd]
        dvelb = AU * [dxbd, dyabd, dzabd]
        return dvelh, dvelb
  
    #General precession from epoch dje to deq.
    deqdat = (dje-dcto-dcbes) / dctrop + dc1900
    prema = premat(deqdat,deq, FK4= True)
    
    dvelh = AU * ( prema.T * [dxhd, dyahd, dzahd] )
    dvelb = AU * ( prema.T * [dxbd, dyabd, dzabd] )
    
    return dvelh, dvelb
#end




###############################################################################
if __name__ == '__main__':
    print airtovac(5167.3216)
    print airtovac(6056.125)
    assert abs(airtovac(6056.125) - 6057.8019)<1e-4
    print airtovac([1234, 5167.3216, 5172.6841, 5183.6042])