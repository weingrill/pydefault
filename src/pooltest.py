'''
Created on Apr 8, 2013

@author: jwe <jweingrill@aip.de>
'''

from multiprocessing import Pool
import numpy

def sqrt(x):
    return numpy.sqrt(x)

if __name__ == '__main__':
    from time import clock, time

    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt

    n = 32
    for k in range(4):
        processes = range(n)
        t = numpy.empty(n)
        for p in processes:
            pool = Pool(p+1)
            start = clock()
            roots = pool.map(sqrt, range(100000))
            dur = clock()-start
            t[p] = dur
            print len(roots)
    
        plt.plot(numpy.array(processes)+1,t)

    plt.show()