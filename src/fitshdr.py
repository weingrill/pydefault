#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Feb 6, 2013

@author: jwe
'''
import pyfits
import sys
hdr = pyfits.getheader(sys.argv[1],1)
print hdr.keys