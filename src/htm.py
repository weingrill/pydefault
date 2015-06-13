#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 12, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import math
import numpy as np

class HTMfunc(object):
    '''
    classdocs
    '''
    verbose = True
    Pi = 3.1415926535897932385E0 
    
    sqrt3    = 1.7320508075688772935
    IDSIZE=64
    IDHIGHBIT = 1L << 63
    IDHIGHBIT2 = 1L << 63
    HTMNAMEMAX = 32
    gEpsilon = 1.0E-15
    HTM_INVALID_ID = 1



    def __init__(self, depth):
        '''
        Constructor


        s0 = {1, 5, 2}; # S0
        S_indexes[0] = s0;
        s1 = {2, 5, 3}; # S1
        S_indexes[1] = s1;
        s2 = {3, 5, 4}; # S2
        S_indexes[2] = s2;
        s3 = {4, 5, 1};# S3
        S_indexes[3] = s3;

        '''
        self.depth = depth
        self.anchor= [[ 0.0,  0.0,  1.0], 
                      [ 1.0,  0.0,  0.0], 
                      [ 0.0,  1.0,  0.0], 
                      [-1.0,  0.0,  0.0], 
                      [ 0.0, -1.0,  0.0], 
                      [ 0.0,  0.0, -1.0]]
        self.S_indexes = [[1, 5, 2], [2, 5, 3], [3, 5, 4], [4, 5, 1]]
        self.N_indexes = [[1, 0, 4], [4, 0, 3], [3, 0, 2], [2, 0, 1]]

    def startpane(v1, v2, v3, xin, yin, zin, name):
        iS2 = 0
        iN1 = 1
        iS1 = 2
        iN2 = 3
        iS3 = 4
        iN0 = 5
        iS0 = 6
        iN3 = 7
        if (xin > 0) and (yin >= 0):
            if zin >= 0: baseindex = iN3 
            else: baseindex = iS0
        
        elif (xin <= 0) and (yin > 0):
            if zin >= 0: baseindex = iN2 
            else: baseindex = iS1
        
        elif (xin < 0) and (yin <= 0):
            if zin >= 0: baseindex = iN1 
            else: baseindex = iS2
        
        elif (xin >= 0) and (yin < 0):
            if zin >= 0: baseindex = iN0 
            else: baseindex = iS3
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
    
        
        baseID = bases[baseindex].ID;
        
        tvec = (double[]) anchor[bases[baseindex].v1];
        v1[0] = tvec[0];
        v1[1] = tvec[1];
        v1[2] = tvec[2];
        
        tvec = (double[])anchor[bases[baseindex].v2];
        v2[0] = tvec[0];
        v2[1] = tvec[1];
        v2[2] = tvec[2];
        
        tvec = (double[])anchor[bases[baseindex].v3];
        v3[0] = tvec[0];
        v3[1] = tvec[1];
        v3[2] = tvec[2];
        
        name.append(bases[baseindex].name);
        return baseID;
    }        

    def m4_midpoint(self, v1, v2):
        w = v1 + v2
        tmp = math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2]*w[2])
        w /= tmp
        return w

    def lookup(self, x, y, z):
        
        
        name= ''
        
        v1 = np.zeros(3)
        v2 = np.zeros(3)
        v0 = np.zeros(3)
        w1 = np.zeros(3)
        w2 = np.zeros(3)
        w0 = np.zeros(3)
        p =  np.zeros(3)
        
        
        p[0] = x;
        p[1] = y;
        p[2] = z;

        # Get the ID of the level0 triangle, and its starting vertices


        startid = self.startpane(v0, v1, v2, x, y, z, name);

        """
        Start searching for the children
        """
        while(self.depth > 0):
            self.depth -= 1
            w2 = self.m4_midpoint(v0, v1, w2)
            w0 = self.m4_midpoint(v1, v2, w0)
            w1 = self.m4_midpoint(v2, v0, w1)
    
            if self.isinside(p, v0, w2, w1):
                name += '0'
                self.copy_vec(v1, w2)
                self.copy_vec(v2, w1)
            elif self.isinside(p, v1, w0, w2):
                name += '1'
                self.copy_vec(v0, v1);
                self.copy_vec(v1, w0);
                self.copy_vec(v2, w2);
            elif self.isinside(p, v2, w1, w0):
                name += '2'
                self.copy_vec(v0, v2);
                self.copy_vec(v1, w1);
                self.copy_vec(v2, w0);
            elif self.isinside(p, w0, w1, w2):
                name += '3'
                self.copy_vec(v0, w0);
                self.copy_vec(v1, w1);
                self.copy_vec(v2, w2);
            else:
                raise(ValueError)
        return name
    
    def radecToVector(self, ra, dec):
        Epsilon = 1.0E-15
        
        vec = np.zeros(3)
        cd = math.cos( (dec * math.pi) / 180.0)

        diff = 90.0 - dec

        if (diff < Epsilon and diff > -Epsilon):
            vec = np.array([1.0, 0.0,1.0])
            return vec

        diff = -90.0 - dec
        
        if (diff < Epsilon and diff > -Epsilon):
            vec = np.array([1.0, 0.0, -1.0])
            return vec;
        
        vec[2] = math.sin((dec* math.pi) / 180.0)
        quadrant = ra / 90.0 # how close is it to an integer?
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
        qint = round(quadrant)
        
        if abs(qint - quadrant) < Epsilon:
            iint = int(qint)
            iint %= 4;
            if (iint < 0): iint += 4

            if iint==0:
                vec[0:2] = [1.0, 0.0]
            elif iint == 1:
                vec[0:2] = [0.0, 1.0]
            elif iint == 2:
                vec[0:2] = [-1.0, 0.0]
            elif iint == 3:
                vec[0:2] = [0.0, -1.0]
            return vec;
        
        vec[0] = math.cos((ra * math.pi) / 180.0) * cd;
        vec[1] = math.sin((ra * math.pi) / 180.0) * cd;
        return vec;

    
#    def lookup(self, ra, dec):
#        v = self.radecToVector(ra,dec)
#        return v
    
    def lookupId(self, ra, dec):
        Pr = math.pi/180.0
        cd = math.cos(dec * Pr);
        x = math.cos(ra * Pr) * cd;
        y = math.sin(ra * Pr) * cd;
        z = math.sin(dec * Pr);
        return np.array([x,y,z])   
    
h = HTMfunc(depth = 10) 
print h.lookup(123.5, -5.123)
print h.lookupId(123.5, -5.123)