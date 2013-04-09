'''
Created on Apr 8, 2013

@author: jwe <jweingrill@aip.de>
'''

def ppsd(t, y):
    from multiprocessing import Pool

    from numpy import linspace, array, cos, sin, dot
    from numpy.linalg import inv

    N = len(y)
    
    def init(argst, argsy):
        global t_global
        global y_global
        
        t_global = argst
        y_global = argsy
        print "init"
    
    def worker(wi):
        global t_global
        global y_global
        
        N = len(t_global)
        print N
        A = array([[cos(wi*ti),sin(wi*ti)] for ti in t_global])
        AT = A.T
        R = dot(AT, A) 
        r = dot(AT, y_global)
        return dot(dot(r.T,inv(R)),r)/N

    w = 2.0*pi*linspace(0.001, 0.5, N)

    pool = Pool(initializer=init, initargs=(t,y))
    p = pool.map(worker, w)
    pool.close() # no more tasks
    pool.join()  # wrap up current tasks
    return p, w/(2.0*pi)

def mpsd(t, y):
    from multiprocessing import Queue, Process

    from numpy import linspace, array, empty, cos, sin, dot, argsort
    from numpy.linalg import inv
    
    def worker(tasks, results):
        # get a task
        task = tasks.get()
        wi = task['wi']
        t = task['t']
        y = task['y']
        N = len(y)
        # do operation
        A = array([[cos(wi*ti),sin(wi*ti)] for ti in t])
        AT = A.T
        R = dot(AT, A) 
        r = dot(AT, y)
        result = dot(dot(r.T,inv(R)),r)/N

        # put the result in results queue
        results.put([wi/(2.0*pi),result])
    

    N = len(y)
    
    w = 2.0*pi*linspace(0.001, 0.5, N)
    f_result = empty(len(w))
    P_result = empty(len(w))

    # create my tasks and results Queues.
    myTasks = Queue()
    myResults = Queue()
    Workers = [Process(target=worker, args=(myTasks, myResults)) for i in range(N)]
    
    for each in Workers:
        each.start()
        
    for i in range(N):
        task = {'wi':w[i], 't':t, 'y':y}
        myTasks.put(task)
        
    while N:
        f_result[N-1],P_result[N-1] = myResults.get()
        N -= 1
    k = argsort(f_result)    
    return [P_result[i] for i in k], [f_result[i] for i in k]   

def psd(t, y):
    """
    Plain Least-Squares Periodogram
    adpted from http://www.sal.ufl.edu/eel6537_2010/LSP.pdf
    """
    from numpy import linspace, array, empty, cos, sin, dot
    from numpy.linalg import inv
    
    N = len(y)
    
    w = 2.0*pi*linspace(0.001,0.5,N)
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
    from numpy import cos, pi
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

    px, f = ppsd(n,y)
    #px, f = psd(n,y)

    
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt
    
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    ax.set_xlabel('n')
    ax.set_ylabel('$y_n$')
    ax.plot(n, y, '.')
    ax.plot(n, y, '--')
    
    ax = fig.add_subplot(2,1,2)
    ax.set_xlabel('f')
    ax.set_ylabel('P($\omega$)')
    ax.plot(f, px)
    
    plt.show()
    pass