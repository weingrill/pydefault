'''
Created on 07.10.2013

@author: jwe
adapted from 
http://math.stackexchange.com/questions/491541/
'''


def interpol(x, y, u):
    p = 0.0
    n = range(len(x))
    for k in n:
        f = 1.0;
        for m in n:
            if m<>k:
                s = np.sin(0.25*(u   -x[m])/np.pi);
                s = s/np.sin(0.25*(x[k]-x[m])/np.pi);
                f *= s;
        p += f*y[k];
    return p

class interp1d(object):
    def __call__(self, x_new):
        
        x_new = np.array(x_new)
        y_new = np.empty(len(x_new))
        i = 0
        for u in x_new:
            p = 0.0
            n = range(len(x))
            for k in n:
                f = 1.0;
                for m in n:
                    if m<>k:
                        s = np.sin(0.25*(u   -x[m])/np.pi);
                        s = s/np.sin(0.25*(x[k]-x[m])/np.pi);
                        f *= s;
                p += f*y[k];
            y_new[i] = p
            i += 1
        return y_new
    
    def __init__(self, x, y):
        import numpy as np
        
        self.x = np.array(x)
        self.y = np.array(y)
    
if __name__ == '__main__':
    import numpy as np
    x = np.array([0,   1,   2,   5,   6,   7])
    y = np.array([0.1, 0.3, 0.4, 0.3, 0.2, 0])
    u = 3.0
    print 'u p(u) = %6.2f %10.2f' % (u, interpol(x, y, u))
    u = 4.0
    print 'u p(u) = %6.2f %10.2f' % (u, interpol(x, y, u))
    
    f = interp1d(x,y)
    print f([3.0,4.0])
    