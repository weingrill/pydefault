#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Sep 18, 2013

@author: jwe
'''
def convert(fitsfile):
    from functions import scaleto
    from PIL import Image
    import numpy as np
    import pyfits
    from os.path import exists

    fileout = os.path.splitext(fitsfile)[0]+'.png'
    if exists(fileout):
        print '%s already exists' % fileout
        return
    
    try:
        hdulist = pyfits.open(fitsfile)
        hdr = hdulist[1].header       
        img = hdulist[1].data
    except IOError:
        exit()
    finally:
        hdulist.close()
    background = hdr['BACKGRND']
    backrms = hdr['BACKRMS']
    warning = ''
    if hdr['MOONDIST']<45.0:
        warning += '%.1fd, ' % hdr['MOONDIST']
    if hdr['MOONILLU']>0.5:
        warning += '%.1f%%, ' % (float(hdr['MOONILLU']) * 100.0)
        
    fimg = np.flipud(img)
    cimg = np.clip(fimg, background/2.0, 65535.0)
    limg = np.log10(cimg)
    simg = scaleto(limg,[0.0,255.0])
    
    a = np.rint(simg)
    #b = a.view('uint8') does not work correct?
    b = a.astype('uint8')
    im = Image.fromarray(b)

    print '%s %8.2f %7.2f %s' % (fitsfile, background, backrms, warning)
    
    im.save(fileout)                    

if __name__ == '__main__':
    import sys, os
    
    if len(sys.argv)>=2:
        for f in sys.argv[1:]:
            convert(f)
    else:
        print 'usage: fits2png infile(s)'