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
    k = 8
    j = np.arange(-k*N, k*N)
    f = 2*j*df/k
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
    k = 16
    j = np.arange(-k*N, k*N)
    f = 2*j*df/k
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

def clean(t, x, gain = 0.5, threshold = 1e-4, maxiter = 1000):
    x -= np.mean(x)
    t -= t[0]
    
    ft, f = dft(t, x)
    ftw, fw = window(t, x)
    aft = abs(ft)
    cleaned = np.zeros(len(ft))
    n2 = len(ft) / 2
    
    for _ in range(maxiter):
        i = np.argmax(aft[n2:])
        ipos = i + n2
        ineg = n2 - i
        amp = gain*aft[ipos]
        sigma = np.std(aft[n2:])
        if amp < sigma or amp < threshold: break
        #print k, i, amp, sigma
        
        cleaned[ipos] += amp/gain
        cleaned[ineg] += amp/gain
        pftw = np.roll(abs(ftw), i)[n2:len(ft)+n2]
        pftw += np.roll(abs(ftw), -i)[n2:len(ft)+n2]
#        plt.plot(f,pftw)
        aft = abs(aft - amp*pftw)
    
    from scipy import signal
    gauss = signal.gaussian(len(cleaned), 2)
    cleaned = signal.convolve(cleaned, gauss, mode='same')
    return f, cleaned+aft, aft

def plot_lightcurve(t, x):
    plt.plot(t, x*1000, 'ko-')
    plt.minorticks_on()
    plt.ylim(max(x*1000)+2, min(x*1000)-2)
    plt.xlim(t[0],t[-1])
    plt.ylabel('mmag')
    
def plot_dft(t, x):
    ft, f = dft(t, x)
    plt.plot(1./f[f>0.0], abs(ft)[f>0.0]/(2.0*np.var(x)), 'k')
    plt.axvline(1, color='r', alpha=0.5)
    plt.axvline(p0, color='b')
    i = np.argmax(ft)
    period = 1./f[i]
    plt.axvline(period, color='g', alpha=0.5)

    
    plt.minorticks_on()
    plt.xlim(0.0,max(t)/3)
    plt.ylabel('S(p)')
    plt.text(0.95, 0.9, 'DFT', va='top', horizontalalignment='right', transform=axis.transAxes)

def plot_window(t, x):
    #window
    ftw, fw = window(t, x)
    plt.plot(1./fw[fw>0.0], abs(ftw)[fw>0.0], 'k')
    plt.xlim(0.0,max(t)/2)
    plt.text(0.95, 0.9, 'window', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)

def plot_lomb(t, x):
    #lomb scargle
    import scipy.signal as signal
    nout = 1000
    flomb = np.linspace(1./60., 1./0.05, nout)
    pgram = signal.lombscargle(t, x, flomb)
    norm = np.sqrt(4*(pgram/nout))/(2.0*np.var(x))
    p = 2.*np.pi/flomb
    plt.plot(p, norm,'k')

    i = np.argmax(norm)
    period = p[i]
    plt.axvline(period, color='g', alpha=0.5)
    plt.axvline(1, color='r', alpha=0.5)
    plt.axvline(p0, color='b', alpha=0.5)

    
    plt.minorticks_on()
    plt.xlim(0.0,max(t)/3)
    plt.ylabel('S(p)')
    plt.text(0.95, 0.9, 'Lomb-Scargle', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)

def plot_clean(t, x):
    f, cleaned, res = clean(t, x, gain=0.9, threshold=2e-3)
    n2 = len(f) /2
    cf = cleaned[n2+1:]/(2.0*np.var(x))
    p = 1./f[n2+1:]
    cf = cf
    p = p
    i = np.argmax(cf)
    period = p[i]
    
    plt.axvline(1, color='r', alpha=0.5)
    plt.axvline(p0, color='b', alpha=0.5)
    plt.axvline(period, color='g', alpha=0.5)
    plt.xlim(0.0,max(t)/3)
    plt.plot(p, cf, 'k')
    plt.minorticks_on()
    plt.ylabel('S(p)')
    plt.text(0.95, 0.9, 'CLEAN', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)
    return period

def plot_phase(t, x, period):
    from functions import phase
    tp, yp = phase(t, x, period)
    plt.scatter(tp,yp*1000, c='k')
    plt.scatter(tp+period,yp*1000, c='k')
    plt.xlim(0.0,period*2)
    plt.ylim(min(yp)*1000,max(yp)*1000)
    plt.xlabel('period = %.2f' % period)
    plt.axvline(period, linestyle='--', color='k')
    plt.minorticks_on()
    plt.ylabel('mmag')
    plt.ylim(plt.ylim()[::-1])
    
def do_shuffle(t, x):
    from random import shuffle
    shuffle(x)
    return t, x
    

def detrend(t, x):
    par = np.polyfit(t, x, 1)
    x -= par[0]*t + par[1]
    return t, x

if __name__ == '__main__':
    from glob import glob
    filenames = glob('/work2/jwe/SOCS/M48/lightcurves.new/*.dat')
    filename = filenames[47]
    filename = '/work2/jwe/SOCS/M48/lightcurves/20140303A-0074-0013#1952.dat'
    t, x, corr = np.loadtxt(filename, unpack=True)
    t -= t[0]
    #x = x[t>30]
    #t = t[t>30]
    #t -= t[0]
    #x -= corr
    x -= np.mean(x)

    t, x = detrend(t, x)
    #t = np.linspace(0,63, 100)
    p0 = 3.24
    #x = 0.3*np.sin(2.*np.pi*t/p0)+0.2
    #x += np.random.normal(scale=0.3, size=x.shape)

    from matplotlib import rcParams
    params = {'backend': 'Agg',
      'axes.labelsize': 8,
      'axes.titlesize': 10,
      'font.size': 8,
      'xtick.labelsize': 8,
      'ytick.labelsize': 8,
      'figure.figsize': [17/2.54, 24/2.54],
      'savefig.dpi' : 300,
      'font.family': 'sans-serif',
      'axes.linewidth' : 0.5,
      #'xtick.major.size' : 2,
      #'ytick.major.size' : 2,
      }
    rcParams.update(params)
    
    plt.figure()
    plt.subplot('511')
    plt.title(filename[32:-4])
    plot_lightcurve(t, x)
    
    axis = plt.subplot('512')
    plot_dft(t, x)
    #sigma = np.std(ft[f>=0])
    
    axis = plt.subplot('513')
    plot_lomb(t, x)

    #lomb scargle
    #nout = 1000
    #flomb = np.linspace(1./60., 1./0.05, nout)
    #freq, pgram = lomb(t, x+14, RMS=1)
    #norm = np.sqrt(4*(pgram/nout))
    #plt.plot(1./freq, pgram,'k')
    #plt.xlim(0.0,max(t)/2)

    
    #from psd import ppsd
    #pp, fp = ppsd(t, x, lower=1./60., upper=1./0.05, num=len(x)*8) 
    #pp = np.sqrt(pp)
    #plt.plot(1./fp, pp,'k')
    
    #plt.xlim(0.0,max(t)/2)
    
    
    #aft = abs(ft)
    
    axis = plt.subplot('514')
    period = plot_clean(t, x)

    plt.subplot('515')
    plot_phase(t, x, period)
    
        
    plt.savefig('/home/jwe/Pictures/Figures/clean.pdf')
    plt.close()
    
    