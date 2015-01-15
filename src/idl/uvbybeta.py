'''
Created on Jan 15, 2015

@author: jwe
'''
from poly import *

def deredd(Eby, by, m1, c1, ub):
    from numpy import array, float64
    Rm1,Rc1, Rub = (-0.33, 0.19, 1.53)
    if type(Eby) is float64: 
        Eby0 =  Eby if Eby>0 else 0.0 
    else:
        Eby0 = array([e if e>0 else 0.0 for e in Eby])
     
    by0 = array(by) - Eby0
    m0 = array(m1) - Rm1*Eby0
    c0 = array(c1) - Rc1*Eby0
    ub0 = array(ub) - Rub*Eby0
    return by0, m0, c0, ub0

def uvbybeta(xby, xm1, xc1, xHbeta, xn, \
             eby_in = None, name = None, verbose=None):
    """
NAME:
      UVBYBETA
PURPOSE:
      Derive dereddened colors, metallicity, and Teff from Stromgren colors.
EXPLANATION:
      Adapted from FORTRAN routine of same name published by T.T. Moon, 
      Communications of University of London Observatory, No. 78. Parameters 
      can either be input interactively (with /PROMPT keyword) or supplied 
      directly.   

CALLING SEQUENCE:
      uvbybeta, /PROMPT               ;Prompt for all parameters
      uvbybeta,by,m1,c1,Hbeta,n        ;Supply inputs, print outputs
      uvbybeta, by, m1, c1, Hbeta, n, Te, Mv, Eby, delm0, radius, 
                      [ TEXTOUT=, Eby_in =, Name =  ]

INPUTS:
      by - Stromgren b-y color, scalar or vector
      m1 - Stromgren line-blanketing parameter, scalar or vector
      c1 - Stromgren Balmer discontinuity parameter, scalar or vector
      Hbeta - H-beta line strength index.  Set  Hbeta to 0 if it is not 
           known, and UVBYBETA will estimate a value based on by, m1,and c1.
           Hbeta is not used for stars in group 8.
      n -  Integer (1-8), scalar or vector,  giving approximate stellar 
           classification

      (1) B0 - A0, classes III - V, 2.59 < Hbeta < 2.88,-0.20 <   c0  < 1.00
      (2) B0 - A0, class   Ia     , 2.52 < Hbeta < 2.59,-0.15 <   c0  < 0.40
      (3) B0 - A0, class   Ib     , 2.56 < Hbeta < 2.61,-0.10 <   c0  < 0.50
      (4) B0 - A0, class   II     , 2.58 < Hbeta < 2.63,-0.10 <   c0  < 0.10
      (5) A0 - A3, classes III - V, 2.87 < Hbeta < 2.93,-0.01 < (b-y)o< 0.06
      (6) A3 - F0, classes III - V, 2.72 < Hbeta < 2.88, 0.05 < (b-y)o< 0.22
      (7) F1 - G2, classes III - V, 2.60 < Hbeta < 2.72, 0.22 < (b-y)o< 0.39
      (8) G2 - M2, classes  IV _ V, 0.20 < m0   < 0.76, 0.39 < (b-y)o< 1.00

OPTIONAL INPUT KEYWORD:
      Eby_in - numeric scalar specifying E(b-y) color to use.   If not
            supplied, then E(b-y) will be estimated from the Stromgren colors
      NAME - scalar or vector string giving name(s) of star(s).  Used only 
              when writing to  disk for identification purposes.
      /PROMPT - if set, then uvbybeta.pro will prompt for Stromgren indicies
               interactively
      TEXTOUT  -  Used to determine output device.  If not present, the
              value of the !TEXTOUT system variable is used (see TEXTOPEN)
              textout=1       Terminal with /MORE (if a tty)
              textout=2       Terminal without /MORE
              textout=3       uvbybeta.prt   (output file)
              textout=4       Laser Printer 
              textout=5       User must open file         
              textout=7       Append to existing uvbybeta.prt file
              textout = filename (default extension of .prt)
     /PRINT - if set, then force display output information to the device 
              specified by !TEXTOUT.    By default, UVBYBETA does not display
              information if output variables are supplied (and TEXTOUT is
              not set). 

OPTIONAL OUTPUTS:
      Te - approximate effective temperature
      MV - absolute visible magnitude
      Eby - Color excess E(b-y)
      delm0 - metallicity index, delta m0, (may not be calculable for early
              B stars).
      radius - Stellar radius (R/R(solar))
EXAMPLE:
      Suppose 5 stars have the following Stromgren parameters

      by = [-0.001 ,0.403, 0.244, 0.216, 0.394 ]
      m1 = [0.105, -0.074, -0.053, 0.167, 0.186 ]
      c1 = [0.647, 0.215, 0.051, 0.785, 0.362] 
      hbeta = [2.75, 2.552, 2.568, 2.743, 0 ]
      nn = [1,2,3,7,8]              ;Processing group number

      Determine stellar parameters and write to a file uvbybeta.prt
      IDL> uvbybeta, by,m1,c1,hbeta, nn, t=3
           ==> E(b-y) = 0.050    0.414   0.283  0.023  -0.025
               Teff =   13060    14030   18420  7250    5760
               M_V =    -0.27    -6.91   -5.94  2.23    3.94
               radius=  2.71     73.51    39.84 2.02    1.53
SYSTEM VARIABLES:
      The non-standard system variables !TEXTOUT and !TEXTUNIT will be  
      automatically defined if they are not already present.   

      DEFSYSV,'!TEXTOUT',1
      DEFSYSV,'!TEXTUNIT',0

NOTES:
      (1) **This procedure underwent a major revision in January 2002
      and the new calling sequence may not be compatible with the old** (NAME
      is now a keyword rather than a parameter.)

      (2) Napiwotzki et al. (1993, A&A, 268, 653) have written a FORTRAN
          program that updates some of the Moon (1985) calibrations.  These
          updates are *not* included in this IDL procedure.
PROCEDURES USED:
      DEREDD, TEXTOPEN, TEXTCLOSE
REVISION HISTORY:                                           
      W. Landsman          IDL coding              February, 1988
      Keyword textout added, J. Isensee, July, 1990
      Made some constants floating point.   W. Landsman    April, 1994
      Converted to IDL V5.0   W. Landsman   September 1997
      Added Eby_in, /PROMPT keywords, make NAME a keyword and not a parameter
                W. Landsman      January 2002
    """
    import numpy as np
    
    Rm1, Rc1, Rub = (-0.33, 0.19, 1.53)
    xby = np.array(xby)
    xm1 = np.array(xm1)
    xc1 = np.array(xc1)
    xHbeta = np.array(xHbeta)
    
    nstar  = xby.size
    xub = xc1 + 2*(xm1 + xby)
    if nstar == 1:
        xflag1 = xHbeta is None
    else:
        xflag1 = [xhb is None for xhb in xHbeta]
    
    do_Eby = eby_in is None
    Te = np.zeros(nstar)
    MV = np.zeros(nstar)
    delm0 = np.zeros(nstar)
    radius = np.zeros(nstar)
    if name is None: 
        name = ['%2d' % (star+1) for star in range(nstar)]
    if not do_Eby: Eby = np.ones(nstar) * eby_in
    else: Eby = np.zeros(nstar)
    for i in range(nstar):
        by = xby[i]
        m1 = xm1[i]
        c1 = xc1[i]
        Hbeta = xHbeta[i]
        n = int(xn[i])
        ub = xub[i]
        flag1 = xflag1[i] 
        flag2 = 0
        warn = ''
        
        if n==1:
            """
            For group 1, beta is a luminosity indicator and c0 is a temperature
            indicator. (u-b) is also a suitable temperature indicator.

            For dereddening a linear relation between the intrinsic (b-y)
            and (u-b) colors is used (Crawford 1978, AJ 83, 48)
            """

            if do_Eby:
                Eby[i] = ( 13.608*by-ub+1.467 ) / (13.608-Rub)
                by0, m0, c0, ub0 = deredd(Eby[i], by, m1, c1, ub)

            """
            If beta is not given it is estimated using a cubic fit to the
            c0-beta relation for luminosity class V given in Crawford (1978).
            """
            if flag1==1: 
                Hbeta = \
                poly(c0, [2.61033, 0.132557, 0.161463, -0.027352] )
            """
            Calculation of the absolute magnitude by applying the calibration
            of Balona & Shobbrock (1974, MNRAS 211, 375)
            """   
            g = np.log10(Hbeta - 2.515) - 1.6*np.log10(c0 +0.322)
            MV[i] = 3.4994 + 7.2026*np.log10(Hbeta - 2.515) -2.3192*g + 2.9375*g**3
            Te[i] = 5040/(0.2917*c0 + 0.2)  

            """
            The ZAMS value of m0 is calculated from a fit to the data of 
            Crawford (1978), modified by Hilditch, Hill & Barnes (1983, MNRAS 204, 241)
            """
            m0zams = poly(c0, [0.07473, 0.109804, -0.139003, 0.0957758] )
            delm0[i] = m0zams - m0
            flag2 = 1
        
        elif n==2:
            if do_Eby:
                """
                For dereddening the linear relations between c0 and (u-b)
                determined from Zhang (1983, AJ 88, 825) is used.
                """
                Eub = ( 1.5*c1 - ub + 0.035) / (1.5/(Rub/Rc1)-1)
                Eby[i] = Eub/Rub
            by0, m0, c0, ub0 = deredd(Eby[i], by, m1, c1, ub)
            if flag1 == 1: Hbeta = 0.037*c0 + 2.542

        elif n==3:
            if do_Eby:
                Eub = (1.36*c1-ub+0.004) / (1.36/(Rub/Rc1)-1)
                Eby[i] = Eub/Rub
            by0, m0, c0, ub0 = deredd(Eby[i], by, m1, c1, ub)
            if flag1: Hbeta = 0.047*c0 +2.578
            
        elif n==4:
            """
            For dereddening the linear relations between c0 and (u-b)
            determined from Zhang (1983, AJ 88, 825) is used.
            """
            if do_Eby:
                Eub = ( 1.32*c1 - ub - 0.056) / ( 1.32 / (Rub/Rc1)-1 )
                Eby[i] = Eub/Rub
                by0, m0, c0, ub0 = deredd(Eby[i], by, m1, c1, ub)
    
            """
            If beta is not given it is derived from a fit of the c0-beta
            relation of Zhang (1983).
            """
            if flag1 == 1: Hbeta = 0.066*c0+2.59
        elif n==5:
            """
            For group 5, the hydrogen Balmer lines are at maximum; hence two
            new parameters, a0 = f{(b-y),(u-b)} and r = f{beta,[c1]} are defined
            in order to calculate absolute magnitude and metallicity.
            """

            if do_Eby:
                m = m1 - Rm1*by
                by0 = 4.2608*m**2 - 0.53921*m - 0.0235
                bycorr = by0
                while ( abs(bycorr - by0) > 0.001):
                    bycorr = by0
                    m0 = m1 - Rm1*(by-bycorr)
                    by0 = 14.0881*m0^2 - 3.36225*m0 + 0.175709
                Eby[i] = by - by0
            by0, m0, c0, ub0 = deredd(Eby[i], by, m1, c1, ub)
            if flag1 == 1: Hbeta = 2.7905 - 0.6105*by + 0.5*m0 + 0.0355*c0
            r = 0.35*(c1-Rc1*by) - (Hbeta-2.565)
            a0 = by0+ 0.18*(ub0-1.36)
            """
            MV is calculated according to Stroemgren (1966, ARA&A 4, 433)
            with corrections by Moon & Dworetsky (1984, Observatory 104, 273)
            """
            MV[i] = 1.5 + 6.0*a0 - 17.0*r
            Te[i] =  5040. /(0.7536 *a0 +0.5282)
            m0zams = -3.95105*by0^2 + 0.86888*by0 + 0.1598
            delm0[i] = m0zams - m0
            
        elif n==6:
            if flag1:
                warn = ' Estimate of Hbeta only valid if star is unreddened'
                Hbeta = 3.06 - 1.221*by - 0.104*c1
    
            m1zams = -2.158*Hbeta**2 +12.26*Hbeta-17.209
            if Hbeta <= 2.74:
                c1zams = 3.0*Hbeta - 7.56
                MVzams = 22.14 - 7*Hbeta

            elif Hbeta > 2.74 and Hbeta <= 2.82 :
                c1zams = 2.0*Hbeta - 4.82
                MVzams = 11.16-3*Hbeta

            else:
                c1zams = 2.0*Hbeta-4.83
                MVzams =-88.4*Hbeta^2+497.2*Hbeta-696.41

            if do_Eby:
                delm1 = m1zams - m1
                delc1 = c1-c1zams
                if delm1 < 0.:
                    by0 = 2.946 - Hbeta - 0.1*delc1 - 0.25*delm1 
                else:
                    by0 = 2.946 - Hbeta - 0.1*delc1
                Eby[i] = by - by0
            by0, m0, c0, ub0 = deredd(Eby[i], by, m1, c1, ub)
            delm0[i] = m1zams - m0
            delc0 = c0 - c1zams
            MV[i] = MVzams -9.0*delc0
            Te[i] = 5040 / (0.771453*by0 + 0.546544)
            
        elif n==7:
            """
            For group 7 c1 is the luminosity indicator for a particular beta,
            while beta {or (b-y)0} indicates temperature.
            Where beta is not available iteration is necessary to evaluate
            a corrected (b-y) from which beta is then estimated.
            """
            if flag1: 
                byinit = by
                m1init = m1
                for _ in range(1000):
                    m1by = 2.5*byinit^2 - 1.32*byinit + 0.345
                    bycorr = byinit + (m1by-m1init) / 2.0
                    if ( abs(bycorr-byinit) <= 0.0001 ): continue
                    byinit = bycorr
                    m1init = m1by
        
                Hbeta = 1.01425*bycorr^2 - 1.32861*bycorr + 2.96618 
    
            """
            m1(ZAMS) and MV(ZAMS) are calculated according to Crawford (1975)
            with corrections suggested by Hilditch, Hill & Barnes (1983,
            MNRAS 204, 241) and Olson (1984, A&AS 57, 443).
            """
            m1zams = poly(Hbeta, [ 46.4167, -34.4538, 6.41701] )
            MVzams = poly(Hbeta, [324.482, -188.748, 11.0494, 5.48012])

            """c1(ZAMS) calculated according to Crawford (1975)"""
            if Hbeta <= 2.65:
                c1zams = 2*Hbeta - 4.91 
            else:
                c1zams = 11.1555*Hbeta**2-56.9164*Hbeta+72.879

            if do_Eby:
                delm1 = m1zams - m1
                delc1 = c1 - c1zams
                dbeta = 2.72 - Hbeta
                by0 = 0.222+1.11*dbeta +2.7*dbeta**2-0.05*delc1-(0.1+3.6*dbeta)*delm1
                Eby[i] = by - by0
            by0, m0, c0, ub0 = deredd(Eby[i], by, m1, c1, ub)
            delm0[i] = m1zams - m0
            delc0 = c0 - c1zams
            f = 9.0 + 20.0*dbeta
            MV[i] = MVzams - f*delc0
            Te[i] = 5040/(0.771453*by0 + 0.546544)

        elif n==8:
            if flag1 == 1: flag1 = 2
            """
            Dereddening is done using color-color relations derived from 
            Olson 1984, A&AS 57, 443)
            """
            if by <= 0.65:
                Eby[i] = (5.8651*by - ub -0.8975) / (5.8651 - Rub)
            elif by > 0.65 and by < 0.79: 
                Eby[i] = (-0.7875*by - c1 +0.6585)/(-0.7875 - Rc1)
                by0 = by - Eby[i]
                if by0 < 0.65:
                    Eby[i] = (5.8651*by - ub -0.8975) / (5.8651-Rub)

            else: 
                Eby[i] = ( 0.5126*by - c1 - 0.3645 ) / (0.5126-Rc1)
                by0 = by - Eby[i]
                if by0 < 0.79: 
                    Eby[i] = (-0.7875*by - c1 + 0.6585) / (-0.7875-Rc1)
                by0  = by - Eby[i]
                if by0 < 0.65: 
                    Eby[i] = ( 5.8651*by - ub - 0.8975) / (5.8651-Rub)

            by0, m0, c0, ub0 = deredd(Eby[i], by, m1, c1, ub)
    

            """
            m1(ZAMS), c1(ZAMS), and MV(ZAMS) are calculated according to Olson (1984)
            """
            m1zams = poly( by0, [7.18436, -49.43695, 122.1875, -122.466, 42.93678]) 
            if by0 < 0.65:
                c1zams = poly(by0, [3.78514, -21.278, 42.7486, -28.7056 ] )
                MVzams =  \
                  poly(by0, [-59.2095, 432.156, -1101.257, 1272.503, -552.48])
            elif (by0 >= 0.65) and (by0 < 0.79):
                c1zams = -0.631821*by0**2 + 0.116031*by0 + 0.33657
                MVzams =   1.37632*by0**2 + 4.97911*by0 + 3.4305
            else:
                c1zams = -0.010028*by0**2 + 0.530426*by0 - 0.37237
                MVzams =  1.18298*by0**2  + 3.92776*by0 + 4.37507
        
            delm0[i] = m1zams - m0
            delc0 =c0 - c1zams
            """
            Teff and MV calibration of Olson (1984)
            """
            if (by0 <= 0.505):
                f = 10. - 80.*(by0-0.38)
                Te[i] = 10**(-0.416*by0+3.924)
            else:
                f = 0.0
                Te[i] = 10**(-0.341*by0+3.869)
        
            MV[i] = MVzams - f*delc0 + 3.2*delm0[i] - 0.07

        else:
            print 'A stellar group of',n,' is not available'
            return
        if n >= 2 and n <= 4:
            """
            c0-beta relation for ZAMS stars according to Crawford (1978,
            AJ 83, 48), modified by Hilditch, Hill & Barnes (1983, MNRAS 204, 241).
            """
            betaza = poly(c0, [2.62745, 0.228638, -0.099623, 0.277363, -0.160402 ] )
            b = betaza - 2.5
            """
            MV(ZAMS) calculated according to Balona & Shobbrock (1984, MNRAS 211, 375)
            """
            MVzams = np.polyval([203.704, -206.98, 77.18, -9.563],b)
            """
            MV is calculated from the d(beta)-d(MV) relation of Zhang (1983)
            """
            dbeta = betaza - Hbeta
            dMV = -121.6*dbeta**2 +61.0*dbeta + 0.08
            MV[i] = MVzams - dMV
            """
            Estimate of Teff by coupling the relations of Boehm-Vitense 
            (1981, ARA&A 19, 295) and Zhang (1983)
            """     
            Te[i] = 5040 / (0.35866*ub0 + 0.27346)
            flag2 = 2
        #endif
        """
        Transformation according to the FV-(b-y)0 relation of Moon 
        (1984, MNRAS 211, 21P)
        """
        if by0 <= 0.335:
            FV = np.polyval([-6.759, 3.731, -1.092, 3.981],by0)
        else: 
            FV = -0.534*by0 + 3.959
        radius[i] = 10.0**(2.*(4.236-0.1*MV[i] - FV))
        if verbose:
            Teff = round(Te[i],-1)
            print '        Star is: ',name[i], \
                    '                Processed in group ', n 
            if len(warn) > 0: print warn
            
            if flag1 == 2: 
                print 'b-y   = %6.3f; m1 = %6.3f; c1 = %6.3f; Hbeta is not used' % \
                    (by,m1,c1)
            elif flag1 == 1: 
                print 'b-y   = %6.3f; m1 = %6.3f; c1 = %6.3f; Hbeta = %5.3f estimated' % \
                    (by, m1,c1, Hbeta) 
            else: 
                print 'b-y   = %6.3f; m1 = %6.3f; c1 = %6.3f; Hbeta = %5.3f' % \
                    (by, m1,c1, Hbeta) 
         
            print '(b-y)0 = %6.3f; m0 = %6.3f: c0 = %6.3f; E(b-y) = %6.3f' % \
                (by0, m0, c0, Eby[i])
        
            print 'Absolute Magnitude (Mv) = %6.2f; Radius (R/R[solar]) = %7.2f' % \
                (MV[i],radius[i])
        
            if flag2 == 2:
                print 'no delta(m0); Effective Temperature (Teff) = %5d K' % Teff 
            else:
                print 'delta(m0) = %6.3f; Effective Temperature (Teff) = %5d K' % \
                    (delm0[i],Teff)
        #end doprint
    #endfor
    return Te,MV,Eby,delm0,radius

if __name__ == '__main__':
    Te,MV,Eby,delm0,radius = uvbybeta(-0.001,0.105,0.647,2.75, 1, verbose=True)
    assert(Eby==0.050)
    by = [-0.001 ,0.403, 0.244, 0.216, 0.394 ]
    m1 = [0.105, -0.074, -0.053, 0.167, 0.186 ]
    c1 = [0.647, 0.215, 0.051, 0.785, 0.362] 
    hbeta = [2.75, 2.552, 2.568, 2.743, None ]
    nn = [1,2,3,7,8]              #Processing group number

    #Determine stellar parameters and write to a file uvbybeta.prt
    Te,MV,Eby,delm0,radius = uvbybeta(by,m1,c1,hbeta, nn, verbose=True)
    print 'E(b-y) = ', '   '.join(['{:.3f}'.format(i) for i in Eby])
    print 'Teff   = ', '   '.join(['{:.0f}'.format(i) for i in Te])
    print 'Mv     = ', '   '.join(['{:.2f}'.format(i) for i in MV])
    print 'radius = ', '   '.join(['{:.2f}'.format(i) for i in radius])
"""
           ==> E(b-y) = 0.050    0.414   0.283  0.023  -0.025
               Teff =   13060    14030   18420  7250    5760
               M_V =    -0.27    -6.91   -5.94  2.23    3.94
               radius=  2.71     73.51    39.84 2.02    1.53

"""