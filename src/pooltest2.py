'''
Created on Apr 9, 2013

@author: jwe <jweingrill@aip.de>
'''

from multiprocessing import Pool
import numpy

def init(argst, argsy):
    global t
    global y
    t = argst
    y = argsy
    print "init"
    
def worker(wi):
    from numpy import array, cos, sin, dot
    from numpy.linalg import inv
    
    global t
    global y
    #y = [0.2,1.0,0.4]

    A = array([[cos(wi*ti),sin(wi*ti)] for ti in t])
    AT = A.T
    R = dot(AT, A) 
    r = dot(AT, y)
    return dot(dot(r.T,inv(R)),r)


def ptest():
    t = [1.0,2.0,3.0]
    y = [0.2,1.0,0.4]
    N = 100

    w = 2.0*numpy.pi*numpy.linspace(0.001, 0.5, N)
    #print worker(w[2])
    pool = Pool(6, initializer=init, initargs=(t,y))
    roots = pool.map(worker, w)
    pool.close() # no more tasks
    pool.join()  # wrap up current tasks

    print roots


if __name__ == '__main__':
    ptest()
