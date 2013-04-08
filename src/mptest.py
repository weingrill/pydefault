'''
Created on Apr 8, 2013

@author: jwe <jweingrill@aip.de>
'''

from multiprocessing import Queue, Process, current_process

def worker(tasks, results):
    # get a task
    t = tasks.get()
    # do operation
    result = t * 2
    # put the result in results queue
    results.put([current_process().name,t,"*",2,"=",result])

if __name__ == '__main__':
    n = 100
    # create my tasks and results Queues.
    myTasks = Queue()
    myResults = Queue()
    Workers = [Process(target=worker, args=(myTasks, myResults)) for i in range(n)]
    
    for each in Workers:
        each.start()
        
    for each in range(n):
        myTasks.put(each)
        
    while n:
        result = myResults.get()
        print "Res: ", result
        n -= 1