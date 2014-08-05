'''
Created on Dec 2, 2013

@author: jwe
'''
def lsqspectrum(t, data, limit=100):
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
    
    amp, freq, phase = calc_ft(t, residual)
    i = np.argmax(amp[:n/2-1])
    fa = amp[i]
    fp = 1./freq[i]
    sigma = np.std(amp[:n/2-1])
    
    func = lambda x, a, p, pha: a*np.cos(2*np.pi/p*x+pha)

    success = 1
    k = 0
    amplitudes = []
    periods = []
    phases = []
    amplitude_errors = []
    period_errors = []
    
    while fa>3.0*sigma and success==1 and k<limit:
        p0 = [fa, fp, phase[i]] # Initial guess for the parameters
        p1, cov = optimize.curve_fit(func, t, residual, p0=p0)

        residual = residual - func(t, p1[0], p1[1], p1[2])
        if p1[0]<0.0:
            p1[0] = -p1[0]
            p1[2]= np.pi + p1[2]
        p1[2] = p1[2] % 2.0*np.pi
        try:
            #calculate error in period
            err = np.sqrt(cov[1,1])
            aerr = np.sqrt(cov[0,0])
        except TypeError:
            err = p1[1]
        if p1[1]>err:
            #print '%7.4f %.4f!' % (p1[1], err)
            amplitudes.append(abs(p1[0]))
            periods.append(p1[1])
            phases.append(p1[2])
            amplitude_errors.append(aerr)
            period_errors.append(err)
            
        else:
            k -= 1
            #print '%7.4f %.4f' % (p1[1], err)
        amp, freq, phase = calc_ft(t, residual)
        i = np.argmax(amp[:n/2-1])
        fa = amp[i]
        if freq[i]>0.0:
            fp = 1./freq[i]
        else:
            fp = 1./freq[i+1]
        sigma = np.std(amp[:n/2-1])
        k += 1
    amplitudes = np.array(amplitudes)
    periods    = np.array(periods)
    phases     = np.array(phases)
    amplitude_errors = np.array(amplitude_errors)
    period_errors = np.array(period_errors)
    residual = np.array(residual
                        )
    return (amplitudes, periods, phases, 
            amplitude_errors, period_errors, 
            residual)

