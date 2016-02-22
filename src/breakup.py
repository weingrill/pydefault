#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jul 30, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

calculate the breakup as a function of mass or B-V respectively
'''
def breakup():
    
    from numpy import sqrt, loadtxt, pi, argmin, arange
    iso = loadtxt('/home/jwe/data/iso_01.000_Gyr.dat')
    #print iso[:,0]
    
    
    mass = iso[:,0]
    radius = iso[:,3]
    bv = iso[:,8]
    logg = iso[:,4]
    
    g = 10**(logg-2) #convert from logg [gcms] to g [kms]
    
    
    minbv = argmin(bv)
    i = arange(minbv)
    
    #G = 6.673e-11
    bv = bv[i]
    radius = radius[i]
    mass = mass[i]
    g = g[i]
    
    
    R_sun = 695660e3
    r = 10**radius * R_sun
    #M_sun = 1.98855e30
    #GM_sun = 1.32712442099e20
    
    #omega = sqrt(GM_sun*mass / r**3)
    
    P = sqrt(4.0*r*pi**2/g)/86400

    for mi,ri,bvi,gi, ppi in zip(mass, radius, bv, logg, P):
        print '%.3f %.3f %.3f %.3f %.3f' % (mi,ri,bvi,gi, ppi)    
    return bv, P

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    bv, P = breakup()
#plt.plot(bv, omega, 'r', label='$\Omega$')
    plt.plot(bv, P,     'g', label='P')
    plt.xlabel('B - V')
    plt.ylabel('P_$breakup$ [days]')
    plt.legend()
    plt.show()
    plt.close()
#print bv
#print P

