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
        import numpy as np
       
        x_new = np.array(x_new)
        y_new = np.empty(len(x_new))
        i = 0
        for u in x_new:
            p = 0.0
            n = range(len(self.x))
            for k in n:
                f = 1.0;
                for m in n:
                    if m<>k:
                        s = np.sin(0.25*(u   -self.x[m])/np.pi);
                        s = s/np.sin(0.25*(self.x[k]-self.x[m])/np.pi);
                        f *= s;
                p += f*self.y[k];
            y_new[i] = p
            i += 1
        return y_new
    
    def __init__(self, x, y):
        import numpy as np
        
        self.x = np.array(x)
        self.y = np.array(y)

def _poolinit(argsx, argsy):
    global px
    global py
    px = argsx
    py = argsy
    
def _poolworker(u):
    from numpy import pi, sin
    
    global px
    global py
    
    p = 0.0
    n = range(len(px))
    fourpi = 4.0*pi
    for k in n:
        f = 1.0;
        for m in n:
            if m<>k:
                f *= sin((u-px[m])/fourpi)/sin((px[k]-px[m])/fourpi);
        p += f*py[k];
    return p

def tinterpol(x, y, x_new):
    from multiprocessing import Pool

    pool = Pool(initializer=_poolinit, initargs=(x, y))
    y_new = pool.map(_poolworker, x_new)
    pool.close() # no more tasks
    pool.join()  # wrap up current tasks

    return y_new

    
if __name__ == '__main__':
    import numpy as np
    x = np.array([0,   1,   2,   5,   6,   7])
    y = np.array([0.1, 0.3, 0.4, 0.3, 0.2, 0])
    u = 3.0
    print 'u p(u) = %6.2f %10.2f' % (u, interpol(x, y, u))
    u = 4.0
    print 'u p(u) = %6.2f %10.2f' % (u, interpol(x, y, u))
    
    f = interp1d(x,y)
    x1 = [3.0,4.0]
    y1 = f(x1)
    x2 = np.linspace(-1.0,8.0, 100)
    y2 =  tinterpol(x, y, x2)
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt
    plt.scatter(x,y)
    plt.scatter(x1,y1,c='g', marker='s')
    plt.plot(x2, y2)
    plt.show()