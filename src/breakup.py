#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jul 30, 2015

@author: Joerg Weingrill <jweingrill@aip.de>

calculate the breakup as a function of mass or B-V respectively
'''
from numpy import sqrt, loadtxt, pi
iso = loadtxt('/home/jwe/data/iso_01.000_Gyr.dat')
#print iso[:,0]


mass = iso[:,0]
radius = iso[:,3]
bv = iso[:,8]
g = iso[:,8]

G = 6.673e-11
bv = bv[mass<2.0]
radius = radius[mass<2.0]
mass = mass[mass<2.0]
g = g[mass<2.0]

R_sun = 695660e3
r = radius * R_sun
M_sun = 1.98855e30
GM_sun = 1.32712442099e20
from matplotlib import pyplot as plt

omega = sqrt(GM_sun*mass / r**3)

#P = sqrt(4.0*r*pi**2/g)/86400
P = (2.0*pi/omega)/86400
plt.plot(bv, P)
plt.show()
#omega = G * M / r**3
print omega