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
    
    hdulist = pyfits.open(fitsfile)
    hdr = hdulist[1].header       
    img = hdulist[1].data
    hdulist.close()
    background = hdr['BACKGRND']
    backrms = hdr['BACKRMS']
    warning = ''
    if hdr['MOONDIST']<45.0:
        warning += '%.1fdeg,' % hdr['MOONDIST']
    if hdr['MOONILLU']>0.5:
        warning += '%.1fper,' % (float(hdr['MOONILLU']) * 100.0)
        
    fimg = np.flipud(img)
    cimg = np.clip(fimg, background/2.0, 65535.0)
    limg = np.log10(cimg)
    simg = scaleto(limg,[0.0,255.0])
    
    a = np.rint(simg)
    #b = a.view('uint8') does not work correct?
    b = a.astype('uint8')
    im = Image.fromarray(b)

    print '%s\t%.2f\t%.2f\t%s' % (fitsfile, background, backrms, warning)
    
    fileout = os.path.splitext(fitsfile)[0]+'.png'
    im.save(fileout)                    

if __name__ == '__main__':
    import pyfits
    import sys, os
    import numpy as np
    
    if len(sys.argv)>=2:
        for f in sys.argv[1:]:
            convert(f)
    else:
        print 'usage: fits2png infile(s)'
                  