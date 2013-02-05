'''
Created on 05.02.2013

@author: jwe
'''
def wcssph2xy(longitude,latitude,x,y,map_type, ctype=None,\
              face=None, pv2 = None,\
              crval=None,crxy=None,longpole=None, latpole= None, \
              north_offset=None, south_offset=None, \
              badindex=None):

# DEFINE ANGLE CONSTANTS
    from math import pi
    from numpy import nan
    pi2 = pi/2.0
    radeg = 57.295779513082323
    map_types=['DEF','AZP','TAN','SIN','STG','ARC','ZPN','ZEA','AIR','CYP',\
            'CAR','MER','CEA','COP','COD','COE','COO','BON','PCO','SFL',\
            'PAR','AIT','MOL','CSC','QSC','TSC','SZP','HPX','HCT']
# check to see that enough parameters (at least 4) were sent
# GENERAL ERROR CHECKING
# find the number of elements in each of the data arrays

    n_long = len( longitude )
    n_lat = len( latitude )
    projection_type = 'DEF'
    if ((map_type >= 0) and (map_type <= 28)):
        projection_type=map_types[map_type]
    else: 
        print 'MAP_TYPE must be >= 0 and < '+\
            len(map_types)+'; it was set to '+\
            str(map_type)

    if projection_type == 'TAN':
        sz_theta = len(theta)
        if sz_theta[0] == 0: 
            x = nan 
        else:
            x = make_array(value = !values.D_NAN, dimen=sz_theta)
        y = x
        g = where(theta > 0)
        if len(g) > 0:
            r_theta = radeg/tan(theta[g])
            x[g] = r_theta*sin(phi[g])
            y[g] = -r_theta*cos(phi[g])
    return

def starast(ra, dec, x, y, cd, righthanded=False, hdr=None, projection=2):    
    '''
    ; NAME:
    ;       STARAST 
; PURPOSE:
;       Compute astrometric solution using positions of 2 or 3 reference stars
; EXPLANATION:
;       Computes an exact astrometric solution using the positions and 
;       coordinates from 2 or 3 reference stars and assuming a tangent 
;       (gnomonic) projection.   If 2 stars are used, then
;       the X and Y plate scales are assumed to be identical, and the
;       axis are assumed to be orthogonal.   Use of three stars will
;       allow a unique determination of each element of the CD matrix.
;
; CALLING SEQUENCE:
;       starast, ra, dec, x, y, cd, [/Righthanded, HDR = h, PROJECTION=]
;
; INPUTS:
;       RA - 2 or 3 element vector containing the Right Ascension in DEGREES
;       DEC- 2 or 3 element vector containing the Declination in DEGREES
;       X -  2 or 3 element vector giving the X position of reference stars
;       Y -  2 or 3 element vector giving the Y position of reference stars
; OUTPUTS:
;       CD - CD (Coordinate Description) matrix (DEGREES/PIXEL) determined 
;               from stellar positions and coordinates.
; OPTIONAL INPUT KEYWORD:
;       /RightHanded - If only 2 stars are supplied, then there is an ambiguity
;               in the orientation of the coordinate system.   By default,
;               STARAST assumes the astronomical standard left-handed system
;               (R.A. increase to the left).   If /Right is set then a 
;               righthanded coordinate is assumed.  This keyword has no effect
;               if 3 star positions are supplied.
;        PROJECTION - Either a 3 letter scalar string giving the projection
;               type (e.g. 'TAN' or 'SIN') or an integer 1 - 25 specifying the
;               projection as given in the WCSSPH2XY procedure.   If not 
;               specified then a tangent projection is computed.
; OPTIONAL INPUT-OUTPUT KEYWORD:
;        HDR - If a FITS header string array is supplied, then an astrometry 
;              solution is added to the header using the CD matrix and star 0
;              as the reference pixel (see example).   Equinox 2000 is assumed.
; EXAMPLE:
;        To use STARAST to add astrometry to a FITS header H;
;
;        IDL> starast,ra,dec,x,y,cd       ;Determine CD matrix
;        IDL> crval = [ra[0],dec[0]]      ;Use Star 0 as reference star
;        IDL> crpix = [x[0],y[0]] +1      ;FITS is offset 1 pixel from IDL
;        IDL> putast,H,cd,crpix,crval     ;Add parameters to header
;
;        This is equivalent to the following command:
;        IDL> STARAST,ra,dec,x,y,hdr=h      
;  
; METHOD:
;       The CD parameters are determined by solving the linear set of equations
;       relating position to local coordinates (l,m)
;
;       For highest accuracy the first star position should be the one closest
;       to the reference pixel.
; REVISION HISTORY:
;       Written, W. Landsman             January 1988
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Added /RightHanded and HDR keywords   W. Landsman   September 2000
;       Write CTYPE values into header   W. Landsman/A. Surkov  December 2002
;       CD matrix was mistakenly transpose in 3 star solution
;       Added projection keyword    W. Landsman   September 2003
;       Test for singular matrix W. Landsman  August 2011 
;-'''
    #from math import pi
    from numpy.linalg import solve
    from numpy import array
#    cdr = pi/180.0 
    map_types=['DEF','AZP','TAN','SIN','STG','ARC','ZPN','ZEA','AIR','CYP',\
            'CAR','MER','CEA','COP','COD','COE','COO','BON','PCO','SFL',\
            'PAR','AIT','MOL','CSC','QSC','TSC']

#    iterate = (len(crpix) == 2) and (len(crval) == 0)
    #Default is tangent proj.
    if type(projection) == str:
        map_type  = map_types.indexof(projection)
        print 'ERROR - supplied projection of ' + projection + ' not recognized'
        map_type = map_type[0]
    else: 
        map_type = projection

    nstar = min( [len(ra), len(dec), len(x), len(y)])
    if (nstar != 2) and (nstar != 3):
        print 'ERROR -  Either 2 or 3 star positions required'
    crval1  = [ ra[0], dec[0] ]
    crpix1  = [ x[0], y[0] ]

    # Convert RA, Dec to Eta, Xi
    eta = 0
    xi = 0
    wcssph2xy( ra[1:], dec[1:], eta, xi, map_type, crval = crval1, latpole = 0.0)
    delx1 = x[1] - crpix1[0] 
    dely1 = y[1] - crpix1[1]     
    
    if nstar == 3:
        delx2 = x[2] - crpix1[0] 
        dely2 = y[2] - crpix1[1]
        b = array([eta[0],xi[0],eta[1],xi[1]])
        a = array( [ [delx1, 0, delx2,    0    ], \
                      [dely1, 0,  dely2,    0  ], \
                      [0. , delx1, 0,    delx2    ], \
                      [0    , dely1   , 0. ,dely2] ] )
    else:
        b = float( [eta[0],xi[0]] )
        if righthanded:
            a = array( [ [delx1,dely1], [-dely1,delx1] ] ) 
        else:
            a = array( [ [delx1,-dely1], [dely1,delx1] ] )

    cd = solve(a,b)        #Solve linear equations
    if nstar == 2:
        if righthanded:
            cd = [ [cd[0],cd[1]],[-cd[1],cd[0]] ] 
        else:
            cd = [ [cd[0],cd[1]],[cd[1],-cd[0]] ] 
    else: 
        cd = cd.reshape(2,2).transpose()


#Add parameters to header
#if N_elements(hdr) GT 0 then begin
#        proj = map_types[map_type]
#        make_astr, astr,CD = cd, crval = crval1, crpix = crpix1+1, $
#                   ctype = ['RA---','DEC--'] + proj
#        putast, hdr, astr, equi=2000.0,cd_type=2
#       
# endif
     
    return

