#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Feb 6, 2013

@author: jwe
'''
import pyfits
import sys
try:
    hdr = pyfits.getheader(sys.argv[1],1)
    for c in hdr.cards: print '%-8s = %s / %s' % c[0:3]
except IndexError:
    hdr = pyfits.getheader(sys.argv[1],0)
    for c in hdr.cards: print '%-8s = %s / %s' % c[0:3]
except IOError:
    print '%s not found' % sys.argv[1]
