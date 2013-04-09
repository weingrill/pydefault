'''
Created on Apr 8, 2013

@author: jwe <jweingrill@aip.de>
'''

import multiprocessing as mp
import Queue

def worker(tasks, results):
    # get a task
    try:
        t = tasks.get(block=False)
    except Queue.Empty:
        pass
    else:
        # do operation
        result = t * 2
        # put the result in results queue
        results.put([mp.current_process().name,t,"*",2,"=",result])

if __name__ == '__main__':
    n = 100
    # create my tasks and results Queues.
    myTasks = mp.Queue(n)
    myResults = mp.Queue(n)
    cores = mp.cpu_count()
    print cores,  " cores"
    Workers = [mp.Process(target=worker, args=(myTasks, myResults)) for i in range(cores)]
    
    for each in Workers:
        each.start()
        
    for each in range(n):
        myTasks.put(each, block=False)
        
    while n:
        if not myResults.empty():
            result = myResults.get(block=False)
            print "Res: ", result
            n -= 1