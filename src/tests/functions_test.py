#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jan 19, 2018

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import functions as fx
import numpy as np
import scipy.signal as sci
from matplotlib import pyplot as plt
x = np.linspace(0,20.0,201)

a = 1.0
b = 20.0
c = 0.5
y1 = fx.gauss(x, a, b, c)
y2 = fx.gaussian(x, a, b, c)
y3 = a*sci.gaussian(len(x), 0.5*(len(x)/np.max(x)))

plt.plot(x, y1, 'r', label='gauss')
plt.plot(x, y2, 'g', label='gaussian')
plt.plot(x, y3, 'b', label='scipy')
plt.legend()
plt.show()