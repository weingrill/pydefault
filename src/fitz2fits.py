#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 30, 2014

@author: jwe

convert a list of fitz files to fits files using imcopy
'''

def convert(fitzfile):
    from subprocess import call
    fitsfile = fitzfile.rstrip('z')+'s'
    call(["imcopy", fitzfile+'[1]', fitsfile])

if __name__ == '__main__':
    import sys
    
    if len(sys.argv)>=2:
        for f in sys.argv[1:]:
            print f
            convert(f)
    else:
        print 'usage: fitz2fits fitzfile(s)'