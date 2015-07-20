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
    k = 3
    j = np.arange(-k*N, k*N)
    f = j*df
    A = np.ones([N, 2*k*N], dtype = complex)
    
    A = A * (-2.j * np.pi * f)
    A = x * np.exp(A.transpose() * t)
    X = np.sum(A, axis=1) / N
    assert(X.shape == f.shape)
    return X, f

def idft(t, X, f=None):
    N = np.shape(t)[0]
    df = 1/(2*(t[-1]-t[0]))
    k = 3
    j = np.arange(-k*N, k*N)
    if f is None: f = j*df
    A = np.ones([N, 2*k*N], dtype = complex)
    
    A = A * (2.j * np.pi * f)
    A = np.exp(A.transpose() * t).transpose() * X.transpose()
    x = np.sum(A, axis=1) /(2*k)
    assert( t.shape == x.shape)
    return t, x

def window(t, x):
    N = np.shape(t)[0]
    df = 1/(2*(t[-1]-t[0]))
    k = 4
    j = np.arange(-k*N, k*N)
    f = j*df
    A = np.ones([N, 2*k*N], dtype = complex)
    
    A = A * (-2.j * np.pi * f)
    A = np.exp(A.transpose() * t)
    X = np.sum(A, axis=1) / N
    print X.shape, f.shape
    return X, f

def clean(t, x, g = 0.2):
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

if __name__ == '__main__':
    #filename = '/work2/jwe/SOCS/M48/lightcurves.new/20140422A-0084-0002#1575.dat'
    #t, x = np.loadtxt(filename, unpack=True)
    t = np.linspace(0,63, 50)
    x = 0.3*np.sin(2.*np.pi*t/0.3237)
    x -= np.mean(x)
    t -= t[0]
    
    plt.figure(figsize=(6,9))

    plt.subplot('411')
    plt.plot(t, x, 'o-')
    
    plt.subplot('412')
    ft, f = dft(t, x)
    plt.plot(f, abs(ft))
    
    plt.subplot('413')
    ftw, fw = window(t, x)
    plt.plot(fw, abs(ftw))
    
    
    #t1 = np.linspace(0,t[-1]-t[0], 1000)
    #x1, _ = idft(t, ft, f=f)
    plt.subplot('414')
    t1,x1 = idft(t, ft, f= f)
    plt.plot(t1, x1, ',-')
    plt.scatter(t, x, edgecolor='none', c='g', alpha=0.5)
    plt.savefig('clean.pdf')
    exit()
    
    i = np.argmax(ft)
    print i
    iw = np.argmax(ftw)
    print iw
    k = ft[i]
    clean = ft - k*np.roll(ftw,i-iw)
    #clean = abs(x1)/abs(x0)
    plt.plot(f,clean ,'r')
    #plt.xlim(0.0,15.0)
    plt.savefig('clean.pdf')
    plt.close()
    exit(0)
    c, f = clean(t-t[0], x)
    plt.plot(f, abs(c))
    #plt.plot(f,A)
    plt.show()
    #print c