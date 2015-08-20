#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jul 13, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import numpy as np
from matplotlib import pyplot as plt

def dft(t, x):
    N = np.shape(t)[0]
    df = 1/(2*(t[-1]-t[0]))
    k = 4
    j = np.arange(-k*N, k*N)
    f = j*df/k
    A = np.ones([N, 2*k*N], dtype = complex)
    
    A = A * (-2.j * np.pi * f)
    A = x * np.exp(A.transpose() * t)
    X = np.sum(A, axis=1) / N
    assert(X.shape == f.shape)
    return X, f

def idft(t, X, f=None):
    N = np.shape(t)[0]
    M = np.shape(f)[0]
    assert X.shape == f.shape
    k = 4.0
    if f is None: 
        k = 3
        df = 1/(2*(t[-1]-t[0]))
        j = np.arange(-k*N, k*N)
        f = j*df
    A = np.ones([N, M], dtype = complex)
    A = A * (2.j * np.pi * f)
    A = np.exp(A.transpose() * t).transpose() * X.transpose()
    x = np.sum(A, axis=1) /(2*k)
    assert( t.shape == x.shape)
    return t, x

def window(t, x):
    N = np.shape(t)[0]
    df = 1/(1*(t[-1]-t[0]))
    k = 8
    j = np.arange(-k*N, k*N)
    f = j*df/k
    A = np.ones([N, 2*k*N], dtype = complex)
    
    A = A * (-2.j * np.pi * f)
    A = np.exp(A.transpose() * t)
    X = np.sum(A, axis=1) / N
    assert( X.shape == f.shape)
    return X, f

def clean_crane(t, x, g = 0.2):
    """
    from Patrick C Crane: Solar Physics 203: 381-408, 2001
    """
    N = np.shape(t)[0]
    x -= np.mean(x)
    dtmin = np.min(abs(t - np.roll(t,1)))
    j = np.arange(-N, N+1)                                                # (A7)
    assert(len(j) == 2*N + 1)
    df = 1/(2*(t[-1]-t[0]))                                              # (A11)
    f = j*df
    print 'dmin = ',dtmin
    print 'N =    ',N
    print 'df =   ',df
    f_c = 1./(2*dtmin)
    fmax = N * df
    print 'fmax = ',fmax
    print 'fc   = ', f_c
    #assert(fmax>=f_c)
    
    def real(f):
        re = np.zeros(2*N + 1)
        for i in range(2*N + 1):
            re[i] = np.sum(x * np.cos(2.0 * np.pi * f[i] * t)) / N
        return re
    
    def imag(f):
        im = np.zeros(2*N + 1)
        for i in range(2*N + 1):
            im[i] = np.sum(x * np.sin(2.0 * np.pi * f[i] * t)) / N
        return im
    
    def win(f): 
        #k = np.arange(-2*N, 2*N+1)
        #fk = k * df                                                      # (A13) 
        w = np.zeros(2*N + 1, dtype = complex)
        for i in range(2*N + 1):
            w[i] = np.sum(np.exp( -2.j * np.pi * f[i] * t)) / N
        return w
    #A = np.sqrt(R**2 + I**2)
    X = real(f) - 1.j*imag(f)
    
    c = np.zeros(2*N + 1, dtype = complex)
    # clean iteration
    
    #1 search maximum in R
    #R = real(f)
    R = X
    #plt.subplot('211')
    #plt.plot(f,X,'r')
    W = win(f)
    #plt.subplot('212')
    #plt.plot(W,'g')
    #plt.show()
    for _ in range(N):
        p = np.argmax(abs(real[f >= 0.0])) + N
        Rp = R[p]
        fp = f[p]
        #print p, fp, Rp, np.max(Rp)
        assert(fp>=0.0)
        
        #plt.plot(f, R, 'r')
        #plt.plot(f, g*Rp*win(f - fp), 'g')
        #plt.axvline(fp, linestyle='--')
        #plt.axhline(Rp, linestyle='--')
        #plt.show()
        
        if fp>0:                                     # (A14a)
            R = R - g*(Rp*win(f - fp) + np.conj(Rp)*win(f + fp))
        elif fp == 0:                                # (A14b)
            R -= g * Rp * win(f)
        
        if fp>0:                                     
            c[f == fp] += g * Rp                     # (A15a)
            c[f == -fp] += g * np.conj(Rp)           # (A15b)
        elif fp == 0.0:
            c[f == 0.0] += g * Rp                      # (A15c)
    
    #def C(f):
        
    return c, f

def clean(t, x, gain = 0.5, threshold = 1e-9):
    x -= np.mean(x)
    t -= t[0]
    ft, f = dft(t, x)
    ftw, fw = window(t, x)
    aft = abs(ft)
    cleaned = np.zeros(len(ft))
    n2 = len(ft) / 2

    for k in range(1000):
        i = np.argmax(aft[n2:])
        ipos = i + n2
        ineg = n2 - i
        amp = gain*aft[ipos]
        sigma = np.std(aft[n2:])
        if amp < sigma or amp < threshold: break
        print k, i, amp, sigma
        
        cleaned[ipos] += amp
        cleaned[ineg] += amp
        pftw = np.roll(abs(ftw), i)[n2:len(ft)+n2]
        pftw += np.roll(abs(ftw), -i)[n2:len(ft)+n2]
#        plt.plot(f,pftw)
        aft = abs(aft - amp*pftw)
    
    from scipy import signal
    gauss = signal.gaussian(len(cleaned), 1.5)
    cleaned = signal.convolve(cleaned, gauss, mode='same')
    return f, cleaned, aft

if __name__ == '__main__':
    from glob import glob
    filenames = glob('/work2/jwe/SOCS/M48/lightcurves.new/*.dat')
    filename = filenames[47]
    t, x = np.loadtxt(filename, unpack=True)

    #t = np.linspace(0,63, 100)
    p0 = 17.2357
    x = 0.3*np.sin(2.*np.pi*t/p0)+0.2
    x += np.random.normal(scale=0.3, size=x.shape)
    #from random import shuffle
    #shuffle(x)
    
    plt.figure(figsize=(6,9))
    plt.subplot('511')
    plt.title(filename[36:-4])
    plt.plot(t, x, 'ko-')
    
    plt.subplot('512')
    ft, f = dft(t, x)
    plt.plot(f, abs(ft), 'k')
    sigma = np.std(ft[f>=0])
    
    plt.subplot('513')
    ftw, fw = window(t, x)
    plt.plot(fw, abs(ftw), 'k')
    print 'ftw', len(ftw)
    
    aft = abs(ft)

    f, cleaned, res = clean(t, x)
    n2 = len(f) /2
    cf = cleaned[n2+1:]/sigma
    p = 1./f[n2+1:]
    cf = cf[p>1.1]
    p = p[p>1.1]
    i = np.argmax(cf)
    if cf[i]>1.0:
        period = p[i]
    else:
        period = 0.0
    
    plt.subplot('514')
    plt.plot(f, res, 'k')
    
    plt.subplot('515')

    plt.axvline(1, color='r')
    plt.axvline(p0, color='b')
    plt.axvline(period, color='g')
    plt.xlim(0.0,max(t)/2)
    plt.plot(p[cf>0.01], cf[cf>0.01], 'k')
    if cf[i]>1.0:
        plt.xlabel('period = %.1f' % period)
        
    plt.savefig('clean.pdf')
    plt.close()
    
    