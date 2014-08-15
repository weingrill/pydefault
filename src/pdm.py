'''
Created on Dec 3, 2013

@author: jwe
'''
def _init(argstime, argsmags, nbins):
    """
    store the lightcurve globally for the workers
    """
    
    global time_global
    global mag_global
    global nbins_global
        
    time_global = argstime
    mag_global = argsmags
    nbins_global = nbins
    
def _worker(period):
    """returns the theta function for the given period"""
    def phase(time, period):
        """
        Calculate phases and ordering.
        """
        from numpy import argsort, floor
        phase = time/period - floor(time/period)
        indi = argsort(phase)
        return phase, indi

    def __setUpEquiBlocks(nbins, phases):
        """
        Set up a sequence of equidistant bins. 
        """
        from numpy import arange, where, logical_and
        blockBegin = (arange(nbins) / float(nbins)).tolist()
        blockEnd   = ((arange(nbins) + 1.0) / float(nbins)).tolist()
        Ns = []
        for i in xrange(nbins):
            indi = where(logical_and(phases >= blockBegin[i], phases < blockEnd[i]))[0]
            Ns.append(len(indi))
    
        i = 0
        while(i < len(Ns)):
            if Ns[i] < 3:
                # The bin with number i does not contain enough data points.
                # Remove the border to make the bins larger...
                iPlus = i + 1
                iMinu = i - 1
                NPlus = 0; NMinu = 0
                if iPlus == len(blockEnd):
                    NPlus = 1e100
                else:
                    NPlus = Ns[iPlus]
                if iMinu == -1:
                    NMinu = 1e100
                else:
                    NMinu = Ns[iMinu]
                if NMinu <= NPlus:
                    # Eliminate the border to bin below
                    blockBegin.pop(i); blockEnd.pop(i-1)
                    Ns[i] = Ns[i-1] + Ns[i]; Ns.pop(i-1)
                    # Reiterate
                    i = 0
                    continue
                else:
                    # Eliminate the border to bin behind
                    blockBegin.pop(i+1); blockEnd.pop(i)
                    Ns[i+1] = Ns[i] + Ns[i+1]; Ns.pop(i)
                    # Reiterate
                    i = 0
                    continue
            i += 1
        return blockBegin, blockEnd
   

    def theta(phase, mag, bbegin, bend):
        """
        Calculate the Theta statistics defined by Stellinger '78.
        """
        from numpy import where, logical_and
        meanMag = mag.mean()
        N = len(mag)
        M = len(bbegin)
        sigmaSqr = ((mag - meanMag)**2).sum() / float(N - 1)

        sSqrUp = 0.0; sSqrDown = 0.0
        for i in xrange(len(bbegin)):
            # Points belonging to a chunk
            indi = where(logical_and(phase >= bbegin[i], phase < bend[i]))[0]
            nj = float(len(indi))
            if nj > 1:
                # Variance of individual block
                sigmaj = ((mag[indi] - mag[indi].mean())**2).sum() / float(nj - 1)
                sSqrUp += (nj - 1)*sigmaj
                sSqrDown += nj
        sSqrDown -= M
    
        sSqr = sSqrUp / sSqrDown
        return sSqr / sigmaSqr

    p, indi = phase(time_global, period)
    bbegin, bend = __setUpEquiBlocks(nbins_global, p[indi])
    return theta(p[indi], mag_global[indi], bbegin, bend)


def pdm(time, mag, minperiod=None, maxperiod=None, delta=None, nbins=None):
    """
    phase dispersion minimization algorithm
    
    minimizes the binned and phased data of a time series given 
    by time and mag. periods are searched within minperiod and maxperiod 
    with the timestep delta.
    
    if nbins is not set the squareroot of datapoints is used for the 
    number of bins.
    
    the periods and their respective theta values are returned
    """
    import numpy as np
    
    if nbins is None:
        nbins = int(np.sqrt(len(time)))
    
    if delta is None:
        delta = np.median(time-np.roll(time, 1))
    
    if minperiod is None:
        minperiod = delta*2.0
    if maxperiod is None:
        maxperiod = (max(time)-min(time))/2
    periods = np.arange(minperiod, maxperiod, delta)
    
    
    from multiprocessing import Pool

    pool = Pool(initializer=_init, initargs=(time, mag, nbins))
    thetas = pool.map(_worker, periods)
    pool.close() # no more tasks
    pool.join()  # wrap up current tasks
    return periods, np.array(thetas)

if __name__ == '__main__':
    
    import numpy as np
    #n = arange(30)
    N = 300
    days = 60.0
    n = np.random.random_sample((N,))*days
    n.sort()
    
    p = 4.0
    amp = 0.1
    sigma = (days/p)/N
    snr = amp/sigma
    noise = np.random.random_sample((N,))* snr
    
    y = amp*np.cos(2*np.pi*n/p)+noise

    periods, theta = pdm(n, y, delta = 30.0/N)
    
    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    ax.set_xlabel('n')
    ax.set_ylabel('$y_n$')
    ax.plot(n, y, '.')
    ax.plot(n, y, '--')
    
    ax = fig.add_subplot(2,1,2)
    ax.set_xlabel('P')
    ax.set_ylabel('$\theta$')
    ax.plot(periods, theta)
    
    ax.axvline(x=periods[np.argmin(theta)], color='r')
    ax.axvline(x=p, linestyle='--', color='g')
    
    plt.show()
