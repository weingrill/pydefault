'''
Created on Apr 8, 2013

@author: jwe <jweingrill@aip.de>
'''

def psd(t, y):
    """
    Plain Least-Squares Periodogram
    adpted from http://www.sal.ufl.edu/eel6537_2010/LSP.pdf
    """
    from numpy import linspace, array, empty, cos, sin, dot
    from numpy.linalg import inv
    
    N = len(y)
    
    w = 2.0*pi*linspace(0.001,0.5,N*10)
    result = empty(len(w))
    i = 0
    for wi in w:
        A = array([[cos(wi*ti),sin(wi*ti)] for ti in t])
        AT = A.T
        R = dot(AT, A) 
        r = dot(AT, y)
        result[i] = dot(dot(r.T,inv(R)),r)/N
        i += 1
    return result,w/(2.0*pi)

if __name__ == '__main__':
    from numpy import arange, cos
    from math import pi
    import numpy as np
    #n = arange(30)
    N = 300
    n = np.random.random_sample((N,))*N
    n.sort()
    f1 = 0.25
    f2 = 0.4
    phi1 = 0.0
    phi2 = 0.784357
    y = 3.0*cos(2*pi*f1*n+phi1) + 4.0*cos(2*pi*f2*n+phi2)
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    ax.set_xlabel('n')
    ax.set_ylabel('$y_n$')
    ax.plot(n, y, '.')
    ax.plot(n, y, '--')
    px, f = psd(n,y)
    
    ax = fig.add_subplot(2,1,2)
    ax.set_xlabel('f')
    ax.set_ylabel('P($\omega$)')
    ax.plot(f, px)
    
    plt.show()
    pass