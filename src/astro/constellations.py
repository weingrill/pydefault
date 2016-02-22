#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jan 19, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''

class Constellations(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self._loadfromfile()
    
    def _loadfromfile(self, filename='/work2/jwe/Projects/Astro/Constellations/eq2000.dat'):
        import numpy as np
        self.data = np.genfromtxt(filename, dtype=(float, float, 'S4','S1'), names=['rah', 'dec', 'const', 'node'], delimiter=' ', autostrip=True) #(11,12,3,1)
        print self.data[:3]
        
    def plot(self, name=None):
        if not name is None:
            pass
        
        if type(name) is str:
            pass
        
if __name__ == '__main__':
    const = Constellations()
    const.plot()