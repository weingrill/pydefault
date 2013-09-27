#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 27, 2013

@author: jwe
'''
import pyfits
import sys
hdr = pyfits.getheader(sys.argv[1])
print hdr[sys.argv[2]]