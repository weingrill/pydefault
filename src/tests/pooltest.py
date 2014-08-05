'''
Created on Apr 8, 2013

@author: jwe <jweingrill@aip.de>
'''

from multiprocessing import Pool
import numpy

c = 1.1

def worker(x):
    return numpy.cos(x)+c

if __name__ == '__main__':
    pool = Pool()
    N = 100
    x = 2.0*numpy.pi*numpy.linspace(0.001, 0.5, N)
    roots = pool.map(worker, x)
    print roots
    
