#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Jul 13, 2015

@author: Joerg Weingrill
@contact: jweingrill@aip.de
@copyright: Copyright 2015, Leibniz Institut für Astrophysik Potsdam (AIP)
@license: MIT
@date: 2015-07-13
@version: 1.0.0
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy import signal


def dft(t, x):
    """
    calculates the discrete Fourier transform of x(t) and returns X(f)
    required by clean function
    """
    N = np.shape(t)[0]
    df = 1 / (2 * (t[-1] - t[0]))
    k = 8
    j = np.arange(-k * N, k * N)
    f = 2 * j * df / k
    A = np.ones([N, 2 * k * N], dtype=complex)

    A = A * (-2.j * np.pi * f)
    A = x * np.exp(A.transpose() * t)
    X = np.sum(A, axis=1) / N
    assert (X.shape == f.shape)
    return X, f


def idft(t, X, f=None):
    """
    calculates the inverse discrete Fourier transform of X(f) and returns x(t)
    """
    N = np.shape(t)[0]
    M = np.shape(f)[0]
    assert X.shape == f.shape
    k = 8
    if f is None:
        k = 8
        df = 1 / (2 * (t[-1] - t[0]))
        j = np.arange(-k * N, k * N)
        f = j * df
    A = np.ones([N, M], dtype=complex)
    A = A * (2.j * np.pi * f)
    A = np.exp(A.transpose() * t).transpose() * X.transpose()
    x = np.sum(A, axis=1) * N / 2  # /(2*k)
    assert (t.shape == x.shape)
    return t, x


def window(t, x):
    """
    calculates the window function of x(t)
    required by clean function
    """
    assert (t.shape == x.shape)
    N = np.shape(t)[0]
    df = 1 / (1 * (t[-1] - t[0]))
    k = 16
    j = np.arange(-k * N, k * N)
    f = 2 * j * df / k
    A = np.ones([N, 2 * k * N], dtype=complex)

    A = A * (-2.j * np.pi * f)
    A = np.exp(A.transpose() * t)
    X = np.sum(A, axis=1) / N
    assert (X.shape == f.shape)
    return X, f


def clean_crane(t, x, g=0.2):
    """
    from Patrick C Crane: Solar Physics 203: 381-408, 2001
    """
    N = np.shape(t)[0]
    x -= np.mean(x)
    dtmin = np.min(abs(t - np.roll(t, 1)))
    j = np.arange(-N, N + 1)  # (A7)
    assert (len(j) == 2 * N + 1)
    df = 1 / (2 * (t[-1] - t[0]))  # (A11)
    f = j * df
    print('dmin = ', dtmin)
    print('N =    ', N)
    print('df =   ', df)
    f_c = 1. / (2 * dtmin)
    fmax = N * df
    print('fmax = ', fmax)
    print('fc   = ', f_c)

    # assert(fmax>=f_c)

    def real(f):
        re = np.zeros(2 * N + 1)
        for i in range(2 * N + 1):
            re[i] = np.sum(x * np.cos(2.0 * np.pi * f[i] * t)) / N
        return re

    def imag(f):
        im = np.zeros(2 * N + 1)
        for i in range(2 * N + 1):
            im[i] = np.sum(x * np.sin(2.0 * np.pi * f[i] * t)) / N
        return im

    def win(f):
        # k = np.arange(-2*N, 2*N+1)
        # fk = k * df                                                      # (A13)
        w = np.zeros(2 * N + 1, dtype=complex)
        for i in range(2 * N + 1):
            w[i] = np.sum(np.exp(-2.j * np.pi * f[i] * t)) / N
        return w

    # A = np.sqrt(R**2 + I**2)
    X = real(f) - 1.j * imag(f)

    c = np.zeros(2 * N + 1, dtype=complex)
    # clean iteration

    # 1 search maximum in R
    # R = real(f)
    R = X
    for _ in range(N):
        p = np.argmax(abs(real[f >= 0.0])) + N
        Rp = R[p]
        fp = f[p]
        assert (fp >= 0.0)

        if fp > 0:  # (A14a)
            R = R - g * (Rp * win(f - fp) + np.conj(Rp) * win(f + fp))
        elif fp == 0:  # (A14b)
            R -= g * Rp * win(f)

        if fp > 0:
            c[f == fp] += g * Rp  # (A15a)
            c[f == -fp] += g * np.conj(Rp)  # (A15b)
        elif fp == 0.0:
            c[f == 0.0] += g * Rp  # (A15c)

    return c, f


def clean(t, x, gain=0.1, threshold=1e-4, maxiter=1000):
    """implementation of the CLEAN algorithm after Roberts et al. (1987)
    :INPUTS:
       t -- (array) independent variable data
       x -- (array) dependent variable data
    
    :OPTIONAL INPUT:
        gain -- gain value for the removal of windowd signal
        threshold -- threshold value, when to abort the iterations
        maxiter -- maximum number of iterations; low gain values need more 
                iterations
    :OUTPUTS:    a tuple of three arrays
        f -- frequencies
        cleaned -- cleaned amplitudes spectrum
        aft -- residual spectrum
    :EXAMPLE:
        t = np.linspace(0.0, 20.0, 100)
        noise = 0.1*np.random.randn(100)
        signal = 7.0 + 0.1*np.sin(t*2.0*np.pi/P) + noise
        window = np.sin(t*2.0*np.pi/1.0)
        i = np.where(window>=0.0)
        t = t[i]
        signal = signal[i]
        f, a, a1 = clean(t, signal, gain=0.5)        
    :NOTES:
        this algorithm can't perform miracles. If the signal is in the order of
        the integrated window amplitude by means snr converges to one, the CLEAN 
        algorithm will fail and pick up spurious signals.
    :REQUIREMENTS:    :doc: `numpy`
    """
    import warnings

    x -= np.mean(x)
    t -= t[0]

    ft, f = dft(t, x)
    ftw, _ = window(t, x)
    aft = abs(ft)
    cleaned = np.zeros(len(ft))
    n2 = len(ft) / 2
    k = 0

    for k in range(maxiter):
        i = np.argmax(aft[n2:])
        ipos = i + n2
        ineg = n2 - i
        amp = gain * aft[ipos]
        assert (aft[ipos] == aft[ineg])
        assert (aft[ipos] == np.max(aft[n2:]))
        sigma = np.std(aft[n2:])
        if amp < gain * sigma or amp < threshold: break

        cleaned[ipos] += amp
        cleaned[ineg] += amp
        pftw = np.roll(abs(ftw), i)[n2:len(ft) + n2]
        pftw += np.roll(abs(ftw), -i)[n2:len(ft) + n2]
        aft = abs(aft - amp * pftw)
    if k == maxiter - 1:
        warnings.warn("maxiter reached!")
    if k == 0:
        warnings.warn("no iterations were performed")
    gauss = signal.gaussian(len(cleaned), 2)
    cleaned = signal.convolve(cleaned, gauss, mode='same')
    return f, cleaned + aft, aft


def plot_lightcurve(t, x):
    plt.plot(t, x * 1000, 'ko-')
    plt.minorticks_on()
    plt.ylim(max(x * 1000) + 2, min(x * 1000) - 2)
    plt.xlim(t[0], t[-1])
    plt.ylabel('mmag')


def plot_dft(t, x, axis):
    ft, f = dft(t, x)
    plt.plot(1. / f[f > 0.0], abs(ft)[f > 0.0] / (2.0 * np.var(x)), 'k')
    plt.axvline(1, color='r', alpha=0.5)
    #    plt.axvline(p0, color='b')
    i = np.argmax(ft)
    period = 1. / f[i]
    plt.axvline(period, color='g', alpha=0.5)

    plt.minorticks_on()
    plt.xlim(0.0, max(t) / 3)
    plt.ylabel('S(p)')
    plt.text(0.95, 0.9, 'DFT', va='top', horizontalalignment='right', transform=axis.transAxes)


def plot_window(t, x, axis):
    # window
    ftw, fw = window(t, x)
    plt.plot(1. / fw[fw > 0.0], abs(ftw)[fw > 0.0], 'k')
    plt.xlim(0.0, max(t) / 2)
    plt.text(0.95, 0.9, 'window', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)


def plot_lomb(t, x, axis):
    # lomb scargle
    import scipy.signal as signal
    nout = 1000
    flomb = np.linspace(1. / 60., 1. / 0.05, nout)
    pgram = signal.lombscargle(t, x, flomb)
    norm = np.sqrt(4 * (pgram / nout)) / (2.0 * np.var(x))
    p = 2. * np.pi / flomb
    plt.plot(p, norm, 'k')

    i = np.argmax(norm)
    period = p[i]
    plt.axvline(period, color='g', alpha=0.5)
    plt.axvline(1, color='r', alpha=0.5)
    #    plt.axvline(p0, color='b', alpha=0.5)

    plt.minorticks_on()
    plt.xlim(0.0, max(t) / 3)
    plt.ylabel('S(p)')
    plt.text(0.95, 0.9, 'Lomb-Scargle', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)


def plot_clean(t, x, axis):
    f, cleaned, _ = clean(t, x, gain=0.9, threshold=2e-3)
    n2 = len(f) / 2
    cf = cleaned[n2 + 1:] / (2.0 * np.var(x))
    p = 1. / f[n2 + 1:]
    cf = cf
    p = p
    i = np.argmax(cf)
    period = p[i]

    plt.axvline(1, color='r', alpha=0.5)
    #    plt.axvline(p0, color='b', alpha=0.5)
    plt.axvline(period, color='g', alpha=0.5)
    plt.xlim(0.0, max(t) / 3)
    plt.plot(p, cf, 'k')
    plt.minorticks_on()
    plt.ylabel('S(p)')
    plt.text(0.95, 0.9, 'CLEAN', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)
    return period


def plot_phase(t, x, period):
    from functions import phase
    tp, yp = phase(t, x, period)
    plt.scatter(tp, yp * 1000, c='k')
    plt.scatter(tp + period, yp * 1000, c='k')
    plt.xlim(0.0, period * 2)
    plt.ylim(min(yp) * 1000, max(yp) * 1000)
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
    x -= par[0] * t + par[1]
    return t, x


def loadfromfile(filename=None):
    if filename is None:
        from glob import glob
        filenames = glob('/work2/jwe/SOCS/M48/lightcurves.new/*.dat')
        filename = filenames[47]
    t, x, _ = np.loadtxt(filename, unpack=True)
    t -= t[0]
    x -= np.mean(x)
    return t, x


def simulate():
    t = np.linspace(0, 63, 100)
    p0 = 3.24
    x = 0.3 * np.sin(2. * np.pi * t / p0) + 0.2
    x += np.random.normal(scale=0.3, size=x.shape)
    return t, x


def acorr(t, x):
    p, f = dft(t, x)
    pc = p.conjugate()
    t1 = np.linspace(t[0], t[-1], 1000)
    ti, xi = idft(t1, p * pc, f=f)
    return ti, xi


def plot_acorr(t, x, axis):
    ta, ac = acorr(t, x)
    print(ta.shape, ac.shape)
    plt.plot(ta, ac, 'k')
    plt.axvline(1, color='r', alpha=0.5)
    #    plt.axvline(period, color='g', alpha=0.5)
    plt.xlim(0.0, max(t) / 3)
    plt.minorticks_on()
    plt.ylabel('S(p)')
    plt.text(0.95, 0.9, 'acorr', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)


def plot_finterpol(t, x):
    p, f = dft(t, x)
    from scipy import signal
    n = 1000
    gaussian_window = signal.gaussian(p.shape[0], std=150)
    p = p * gaussian_window
    t1 = np.linspace(t[0], t[-1], n)
    ti, xi = idft(t1, p, f=f)

    print(ti.shape, xi.shape)
    plt.scatter(t, x * 1000, alpha=0.5)
    plt.plot(ti, xi, 'k')
    plt.xlim(ti[0], ti[-1])
    plt.minorticks_on()
    plt.ylabel('mag')
    # plt.text(0.95, 0.9, 'acorr', verticalalignment='top', horizontalalignment='right', transform=axis.transAxes)


def _test():
    filename = '/work2/jwe/SOCS/M48/lightcurves/20140303A-0074-0013#1952.dat'
    t, x = loadfromfile(filename)
    t, x = detrend(t, x)

    from matplotlib import rcParams
    params = {'backend': 'Agg',
              'axes.labelsize': 8,
              'axes.titlesize': 10,
              'font.size': 8,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8,
              'figure.figsize': [17 / 2.54, 24 / 2.54],
              'savefig.dpi': 300,
              'font.family': 'sans-serif',
              'axes.linewidth': 0.5,
              # 'xtick.major.size' : 2,
              # 'ytick.major.size' : 2,
              }
    rcParams.update(params)

    plt.figure()
    plt.subplot('511')
    plt.title(filename[32:-4])
    plot_lightcurve(t, x)

    axis2 = plt.subplot('512')
    plot_dft(t, x, axis2)
    # sigma = np.std(ft[f>=0])

    plt.subplot('513')
    # plot_lomb(t, x)
    # plot_acorr(t, x)
    plot_finterpol(t, x)

    axis4 = plt.subplot('514')
    period = plot_clean(t, x, axis4)

    plt.subplot('515')
    plot_phase(t, x, period)

    plt.savefig('/home/jwe/Pictures/Figures/clean.pdf')
    plt.close()


if __name__ == '__main__':
    _test()
