#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 27, 2013

@author: Joerg Weingrill <jweingrill@aip.de>
'''

if __name__ == '__main__':
    import argparse
    import astropy.io.fits as pyfits
    
    parser = argparse.ArgumentParser(description='reads a specific key of a FITS file')
    parser.add_argument('-key', help='FITS header key')
    parser.add_argument('filename', help='name of FITS file')

    args = parser.parse_args()

    hdr = pyfits.getheader(args.filename)
    print hdr[args.key]