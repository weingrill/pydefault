#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 12, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

Revision History
================

$Log: HTMfunc.java,v $
Revision 1.5  2003/02/19 15:46:11  womullan
Updated comments mainly to make java docs better

Revision 1.4  2003/02/14 21:46:22  womullan
*** empty log message ***

Revision 1.3  2003/02/14 16:46:00  womullan
Newly organised classes in the new packages

Revision 1.2  2003/02/11 18:40:27  womullan
moving to new packages and names

Revision 1.1  2003/02/07 23:02:51  womullan
Major tidy up of HTM to HTMfunc addition of exceptions and error handler

Revision 1.2  2003/02/07 19:36:03  womullan
  new HTM for cc methods

Revision 1.1  2003/02/05 23:00:07  womullan
Added new timimig routines and HTM clas to do cc stuff from the c++ package
'''
import math
import numpy as np

class HTMfunc(object):
    '''
    classdocs
    '''

    def __init__(self, depth):
        '''
        Constructor
        '''
        self.depth = depth
        self.name = ''

    def startpane(self, xin, yin, zin):
        """
        where to start in the HTM quadtree when looking for a vectors home
        xin yin zin are thin input vector
        v1 v2 v3 and name are loaded with the initial triangle points
        and the name of the triangle
        """
        if (xin > 0) and (yin >= 0):
            baseindex = 'N3' if zin >= 0 else 'S0'
        
        elif (xin <= 0) and (yin > 0):
            baseindex = 'N2' if zin >= 0 else 'S1'
        
        elif (xin < 0) and (yin <= 0):
            baseindex = 'N1' if zin >= 0 else 'S2'
        
        elif (xin >= 0) and (yin < 0):
            baseindex = 'N0' if zin >= 0 else 'S3'
        else:
            raise(ValueError)

        bases= {}
        bases["S2"] = [10, 3, 5, 4]
        bases["N1"] = [13, 4, 0, 3]
        bases["S1"] = [9, 2, 5, 3]
        bases["N2"] = [14, 3, 0, 2]
        bases["S3"] = [11, 4,5,1]
        bases["N0"] = [12, 1, 0, 4]
        bases["S0"] = [8, 1, 5, 2]
        bases["N3"] = [15, 2, 0, 1]

        anchor= [np.array([ 0.0,  0.0,  1.0]), 
                 np.array([ 1.0,  0.0,  0.0]), 
                 np.array([ 0.0,  1.0,  0.0]), 
                 np.array([-1.0,  0.0,  0.0]), 
                 np.array([ 0.0, -1.0,  0.0]), 
                 np.array([ 0.0,  0.0, -1.0])]
        
        self.baseID = bases[baseindex][0]
        tvec = anchor[bases[baseindex][1]]
        v1 = np.array(tvec)
        
        tvec = anchor[bases[baseindex][2]]
        v2 = np.array(tvec)
        
        tvec = anchor[bases[baseindex][3]]
        v3 = np.array(tvec)
        
        self.baseindex = baseindex
        self.name = baseindex
        return v1, v2, v3

    def m4_midpoint(self, v1, v2):
        """
        calculate midpoint fo vectors v1 and v2 the answer is put in
        the provided w vector.
        """
        w = v1 + v2
        tmp = np.sqrt(np.sum(w**2))
        w /= tmp
        return w

    def _lookup(self, x, y, z):
        """
        for a given vector x,y,z return the HTM ID to the given depth
        """
        
        w1 = np.zeros(3)
        w2 = np.zeros(3)
        w0 = np.zeros(3)
        
        p = [x, y, z]
        
        # Get the ID of the level0 triangle, and its starting vertices


        v0,v1,v2 = self.startpane(x, y, z)

        """
        Start searching for the children
        """
        depth = self.depth
        while(depth > 0):
            depth -= 1
            w2 = self.m4_midpoint(v0, v1)
            w0 = self.m4_midpoint(v1, v2)
            w1 = self.m4_midpoint(v2, v0)
    
            if self.isinside(p, v0, w2, w1):
                self.name += '0'
                v1 = w2.copy()
                v2 = w1.copy()
            elif self.isinside(p, v1, w0, w2):
                self.name += '1'
                v0 = v1.copy()
                v1 = w0.copy()
                v2 = w2.copy()
            elif self.isinside(p, v2, w1, w0):
                self.name += '2'
                v0 = v2.copy()
                v1 = w1.copy()
                v2 = w0.copy()
            elif self.isinside(p, w0, w1, w2):
                self.name += '3'
                v0 = w0.copy()
                v1 = w1.copy()
                v2 = w2.copy()
            else:
                raise(ValueError)
        return self.name

    def _lookupId(self, x, y, z):
        """
        looks up the name of a vector x,y,z and converts the name to an id
        """
        self.name = self._lookup(x, y, z)
        return self.nameToId()

    def lookup(self, ra, dec):
        """
        For given ra and dec lookup the HTMID to given depth
        HTM works in vectors so this basically converts ra dec to a vector and
        calls lookup for the vector.
        """
        x,y,z = self.radecToVector(ra,dec)
        return self._lookup(x, y, z)
    
    def radecToVector(self, ra, dec):
        """convert ra dec to a vector"""
        Epsilon = 1.0E-7
        quadrant = (ra/90.0) # how close is it to an integer?
        ra *= math.pi / 180.0
        dec *= math.pi / 180.0
        vec = np.zeros(3)
        
        if abs(90.0 - dec) < Epsilon:
            return np.array([1.0, 0.0,1.0])

        if abs(-90.0 - dec) < Epsilon:
            return np.array([1.0, 0.0, -1.0])
        
        vec[2] = math.sin(dec)
        
        """
        if quadrant is (almost) an integer, force x, y to particular
        values of quad:
        quad,   (x,y)
        0       (1,0)
        1,      (0,1)
        2,      (-1,0)
        3,      (0,-1)
        q>3, make q = q mod 4, and reduce to above
        q mod 4 should be 0.
        """
        vec[0] = math.cos(dec) * math.cos(ra)
        vec[1] = math.cos(dec) * math.sin(ra)
        return vec
        qint = round(quadrant)
        
        if abs(qint - quadrant) < Epsilon:
            #print 'q! %d %f' % (qint, quadrant),
            iint = int(qint)
            iint %= 4
            if (iint < 0): iint += 4

            if iint==0:
                vec[0:2] = [1.0, 0.0]
            elif iint == 1:
                vec[0:2] = [0.0, 1.0]
            elif iint == 2:
                vec[0:2] = [-1.0, 0.0]
            elif iint == 3:
                vec[0:2] = [0.0, -1.0]
            return vec
        
        vec[0] = math.cos(dec) * math.cos(ra)
        vec[1] = math.cos(dec) * math.sin(ra)
        return vec
    
    def lookupId(self, ra, dec):
        """
        same as lookup but converts the name to an id
        """
        ra *= math.pi/180.0
        dec *= math.pi/180.0
        x = math.cos(dec) * math.cos(ra)
        y = math.cos(dec) * math.sin(ra)
        z = math.sin(dec)
        return self._lookupId(x, y, z) 
    
    def isinside(self, p, v1, v2, v3): # p need not be normalized!!!
        """
        for a given vector p is it contained in the triangle whose corners are
        given by the vectors v1, v2, v3.
        """
        gEpsilon = 1.0E-15
        
        crossp = np.cross(v1, v2)
        if (np.dot(p, crossp) < -gEpsilon):
            return False

        crossp = np.cross(v2, v3)
        if (np.dot(p, crossp) < -gEpsilon):
            return False

        crossp = np.cross(v3, v1)
        if (np.dot(p, crossp) < -gEpsilon):
            return False

        return True
    
    def nameToId(self):
        """
        * for a given name i.e. N301022 convert it to its 64bit htmId
        *  effectively this walks trough the string in reverse order and
        * sets the bits of a long number according to the charector in the
        * String.
    
           The  name has always the same structure, it begins with
           an N or S, indicating north or south cap and then numbers 0-3 follow
           indicating which child to descend into. So for a depth-5-index we have
           strings like
                     N012023  S000222  N102302  etc
    
           Each of the numbers correspond to 2 bits of code (00 01 10 11) in the
           uint64. The first two bits are 10 for S and 11 for N. For example
    
                     N 0 1 2 0 2 3
                     11000110001011  =  12683 (dec)
    
           The leading bits are always 0
    
        """
        HTMNAMEMAX = 32

        if self.name == None or len(self.name) == 0: # null pointer-name
            raise(ValueError)
        if not self.name[0] in ['N','S']:  # invalid name
            raise(ValueError)

        siz =len(self.name)       # determine string length
        # at least size-2 required, don't exceed max
        if siz < 2:
            raise(ValueError)
        if siz > HTMNAMEMAX:
            raise(ValueError)
        
        out = 0
        for i in range(siz-1, 0, -1):
            # set bits starting from the end
            if self.name[i] > '3' or self.name[i] < '0': 
                # invalid name
                raise(ValueError)

            out +=  int(self.name[i]) << 2*(siz - i -1)

        i = 2 # set first pair of bits, first bit always set
        if self.name[0] == 'N': 
            i += 1 # for north set second bit too
        last = i << (2*siz - 2)
        out += last
        return out
    
    def idLevel(self, htmid):
        """
        idLevel is a trusting method (assumes that the id is well formed and
        valid) that returns the level of the trixel represented by the given
        64 bit htmId.
        This determines the index of the first bit set by continually left
        shifting the htmid 2 bits at a time until a set bit is found in the
        highest position.
        2 bits are used as each number in an HTM Name like N012 is represented
        by two bits.
        Also the number of bits used dived by two is the level of the htmId +2
        we have +2 because of the N0,S1 etc  prefixes.
        So an id with bit 64 set is level 30 (64/2 -2)
        """
        IDSIZE = 64
        IDHIGHBIT = 1L << 63
        
        # determine index of first set bit
        for i in range(0, IDSIZE, 2):
            if  ((htmid << i) & IDHIGHBIT )> 0:
                break
        size=(IDSIZE-i) >> 1
        """
        Size is the length of the string representing the name of the
           trixel, the level is size - 2
        """
        return size-2
    
    def idToName(self, htmid):
        """
        Walk the bits of the id and convert it to a string like N012.
        sucessivley look at each pair of bits and convert it to 0,1,2 or 3.
        Finally use the highest bit to set N or S in the string.
    
        The first bit pair gives N (11) or S (10).
        The subsequent bit-pairs give the numbers 0-3: (00, 01, 10, 11).

        Example: In the first level there are 8 nodes, and we get
        the names from numbers 8 through 15 giving S0,S1,S2,S3,N0,N1,N2,N3.

        The order is always ascending starting from S0000.. to N3333...
        """
        IDSIZE = 64
        IDHIGHBIT = 1L << 63
        IDHIGHBIT2 = 1L << 63
    

        # determine index of first set bit
        for i in range(0,IDSIZE,2):
            x8 = ((htmid << i) & IDHIGHBIT)
            x4 = ((htmid << i) & IDHIGHBIT2)
            if ( x8 != 0 ): 
                break
            if (x4 != 0): # invalid id
                raise(ValueError)
        
        if (htmid == 0): raise(ValueError)
        size=(IDSIZE-i) >> 1

        name = ''
        #fill characters starting with the last one
        for i in range(size-1):
            c =  '%d' % ((htmid >> i*2) & 3)
            name = c + name


        # put in first character
        if ((htmid >> (size*2-2)) & 1) > 0:
            name = 'N' + name
        else:
            name = 'S' + name

        return name
    
    def nameToTriangle(self, name):
        """
        Like name2Id but instaed of returning the htmId
        return v0,v1,v2 vetors representin the corners of the trinagle
        @return Object[] which is holding 3 double[] for v0 v1 v2
        """
        S_indexes = [np.array([1, 5, 2]), 
                     np.array([2, 5, 3]), 
                     np.array([3, 5, 4]), 
                     np.array([4, 5, 1])]
        N_indexes = [np.array([1, 0, 4]), 
                     np.array([4, 0, 3]), 
                     np.array([3, 0, 2]), 
                     np.array([2, 0, 1])]
       
        # Get the top level hemi-demi-semi space
        
        anchor_offsets= np.zeros(3)
        k = int(name[1])

        if name[0] == 'S': 
            anchor_offsets = S_indexes[k]
        else:
            anchor_offsets = N_indexes[k]

        anchor= [np.array([ 0.0,  0.0,  1.0]), 
                 np.array([ 1.0,  0.0,  0.0]), 
                 np.array([ 0.0,  1.0,  0.0]), 
                 np.array([-1.0,  0.0,  0.0]), 
                 np.array([ 0.0, -1.0,  0.0]), 
                 np.array([ 0.0,  0.0, -1.0])]

        v0 = anchor[anchor_offsets[0]].copy()
        v1 = anchor[anchor_offsets[1]].copy()
        v2 = anchor[anchor_offsets[2]].copy()
        offset = 2
        while offset < len(name):
            s = name[offset]
            w2 = self.m4_midpoint(v0, v1)
            w0 = self.m4_midpoint(v1, v2)
            w1 = self.m4_midpoint(v2, v0)
            if s == '0':
                #v0 = v0 
                v1 = w2.copy()
                v2 = w1.copy()
            elif s == '1':
                v0 = v1.copy()
                v1 = w0.copy()
                v2 = w2.copy()
            elif s == '2':
                v0 = v2.copy()
                v1 = w1.copy()
                v2 = w0.copy()
            elif s == '3':
                v0 = w0.copy()
                v1 = w1.copy()
                v2 = w2.copy()
            offset += 1
        return [v0,v1,v2]

    def _idToPoint(self, name):
        """
        for a given ID get back the approximate center of the triangle.
        This may be used as an inverse of lookup however bear in mind many points
        may fall in a triangle while only the center point will be returned from
        this function.
        The name is used to get the triange corners and these are then averaged
        to get a center point.
        The resultant vector is not normalized.
        """
        tri = self.nameToTriangle(name)
        center = np.sum(tri, axis=0)
        csum = np.sqrt(np.sum(center**2))
        center /= csum
        
        return center # I don't want it nomralized or radec to be set,

    def idToPoint(self, htmid):
        """
        gets the name from the id and calls idToPoint with it.
        """
        name = self.idToName(htmid)
        return self._idToPoint(name)

    def _distance(self, v1, v2):
        """
        return the angular distance between two vectors
        ACOS (V1 . V2)
        """
        return math.acos(np.dot(v1,v2))

    def distance(self, htmId1, htmId2):
        """
        return the angular distance between two htmids
        gets the vectors of the mid points and uses thoose to compute distance
        """
        if type(htmId1) is int:               
            v1 = self.idToPoint(htmId1)
            v2 = self.idToPoint(htmId2)
            return self._distance(v1,v2)
        
        elif type(htmId1) is str:
            v1 = self.idToPoint(htmId1)
            v2 = self.idToPoint(htmId2)
            return self._distance(v1,v2)
        
        else:
            raise(TypeError)
    
    def Vectortoradec(self, xvec):
        """
        convert the cartesian coodinates x 
        to spherical coordinates lamb, beta and r
        result in degrees
        taken from Montenbruck p.2
        """
        x, y, z = xvec
        r = np.sqrt(x**2 + y**2 + z**2)
        rho = np.sqrt(x**2 + y**2)
        
        if rho == 0.0:
            if z > 0.0:  beta =  90.0
            if z == 0.0: beta =   0.0
            if z < 0.0:  beta = -90.0
        else:
            beta = math.atan2(z, rho) * 180.0 / math.pi
        
        phi = 2.0 * math.atan2(y, abs(x)+rho) * 180.0 / math.pi
        
        if x == 0.0 and y == 0.0: lamb = 0.0
        if x >= 0.0 and y >= 0.0: lamb = phi
        if x >= 0.0 and y < 0.0:  lamb = 360.0 + phi
        if x < 0.0:               lamb = 180.0 - phi
        
        return lamb, beta, r

if __name__ == '__main__':
    h = HTMfunc(depth = 12) 
    
    ra0 = 45.0; dec0 = 45.0
    htmid = h.lookupId(ra0, dec0)
    print 'baseindex:', h.baseindex
    print 'baseID:', h.baseID
    print 'basename:', h.name
    print 'id = %d' % htmid
    x = h.idToPoint(htmid)
    #x = h._idToPoint(h.name)
    ra, dec, r = h.Vectortoradec(x) 
    #x = [0.707106781187,-0.707106781187, 0.0]
    print '(%5.1f, %5.1f) | (%6.3f, %6.3f, %6.3f) | (%5.1f, %5.1f)' % \
        (ra0, dec0, x[0], x[1], x[2], ra,dec)
    htmid = 12683
    print h.idToName(htmid)
    #N012023
    
    for dec0 in [-90.0,-45.0,0.0,45.0,90.0]:
        for ra0 in range(0, 360, 1):
            x1 = h.radecToVector(float(ra0), dec0)
            ra, dec, r = h.Vectortoradec(x1) 
            if abs(ra-ra0)>1e-7 or abs(dec-dec0)>1e-7:
                print 'ra0, dec0 = %5.1f, %5.1f' % (ra0, dec0),
                print 'x1 = %6.3f, %6.3f, %6.3f' % (x1[0], x1[1], x1[2]),
                print 'ra, dec = %5.1f, %5.1f' % (ra, dec)
            
    
   
    #import esutil
    #h1 = esutil.htm.HTM(8)
    #htmid1 = h1.lookup_id(ra0, dec0)[0] 
    #print htmid1
    
    