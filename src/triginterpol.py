'''
Created on 07.10.2013

@author: jwe
adapted from 
http://math.stackexchange.com/questions/491541/trigonometric-interpolation
'''

import numpy as np
x = np.array([0,   1,   2,   5,   6,   7])
y = np.array([0.1, 0.3, 0.4, 0.3, 0.2, 0])

def interpol(u):
    p = 0.0
    for k in range(len(x)):
        f = 1.0;
        for m in range(len(x)):
            if m<>k:
                s = np.sin(0.25*(u   -x[m])/np.pi);
                s = s/np.sin(0.25*(x[k]-x[m])/np.pi);
                f *= s;
        p += f*y[k];
    return p

if __name__ == '__main__':
    u = 3.0
    print 'u p(u) = %6.2f %10.2f' % (u, interpol(u))
    u = 4.0
    print 'u p(u) = %6.2f %10.2f' % (u, interpol(u))
    
    