'''
Created on Dec 2, 2013

@author: jwe
'''
def lsqspectrum(t, data, limit=100, residual=None):
    """
    compute the least-sqares-spectrum of a given dataset
    """
    import numpy as np
    from math import atan2
    from scipy import optimize  
    
    n = len(data)
    def calc_ft(t, y):
        ft = np.fft.rfft(y, n)
        amp = 2.*abs(ft)/n
        freq = np.fft.fftfreq(n, d=abs(t[1]-t[0]))
        phase = [atan2(c.real, c.imag) for c in ft]
        return (amp, freq, phase)
    
    residual = data
    
    amp,freq,phase = calc_ft(t, residual)
    i = np.argmax(amp[:n/2-1])
    fa = amp[i]
    fp = 1./freq[i]
    sigma = np.std(amp[:n/2-1])
    
    fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/p[1]*x+p[2])
    errfunc = lambda p, x, y: y - fitfunc(p, x) # Distance to the target function

    success = 1
    k = 0
    amplitudes = []
    periods = []
    phases = []
    
    while fa>3.0*sigma and success==1 and k<limit:
        p0 = [fa, fp, phase[i]] # Initial guess for the parameters
        
        p1, cov, _, _, success = optimize.leastsq(errfunc, p0[:], args=(t, residual), full_output=1)
        residual = residual - fitfunc(p1, t)
        if p1[0]<0.0:
            p1[0] = -p1[0]
            p1[2]= np.pi + p1[2]
        p1[2] = p1[2] % 2.0*np.pi
        try:
            err = np.sqrt(cov[1,1])
        except TypeError:
            err = p1[1]
        #print '%7.4f %.4f %.1f' % (p1[1], p1[0], fa/sigma)
        if success and p1[1]>err:
            print '%7.4f %.4f!' % (p1[1], err)
            amplitudes.append(abs(p1[0]))
            periods.append(p1[1])
            phases.append(p1[2])
        else:
            k -= 1
            print '%7.4f %.4f' % (p1[1], err)
        amp, freq, phase = calc_ft(t, residual)
        i = np.argmax(amp[:n/2-1])
        fa = amp[i]
        if freq[i]>0.0:
            fp = 1./freq[i]
        else:
            fp = 1./freq[i+1]
        sigma = np.std(amp[:n/2-1])
        k += 1
    return (amplitudes, periods, phases, residual)

    
    
    
