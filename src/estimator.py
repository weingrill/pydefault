'''
Created on Jan 18, 2013

@author: jwe
'''

from time import clock

class Estimator(object):
    '''
    Estimator object to provide an estimate and a percentage of a process
    '''


    def __init__(self, firstitem, lastitem):
        '''
        Constructor
        '''
        self.firstitem = firstitem
        self.lastitem = lastitem
        self.s0 = clock()
    
    def estimate(self, current, comment=''):
        s1 = clock()
        if current < self.firstitem:
            self.firstitem = current  
        totalitems = self.lastitem - self.firstitem
        percent = 100.*(current-self.firstitem)/totalitems
        perfile = (s1-self.s0)/(current-self.firstitem+1)
        togo = (totalitems - current + self.firstitem)*perfile
        hours = int(togo // 3600)
        mins =  int((togo % 3600) // 60)
        secs = (togo % 3600) % 60
        togostr = '%02d:%02d:%04.1f' % ( hours, mins, secs)
        print "%.2f%%, %.3fs/item, %s: %s" % (percent, perfile, togostr, comment)
