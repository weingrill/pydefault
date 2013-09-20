'''
Created on Sep 18, 2013

@author: jwe
'''
from functions import scaleto
from PIL import Image
import numpy as np
import pyfits
f = '/work2/jwe/stella/wifsip/20130918/science20130917A-0051EXP0009.fitz'
hdulist = pyfits.open(f)
hdr = hdulist[1].header       
img1 = hdulist[1].data
logim = np.log10(img1)
scalim = scaleto(logim,[0,255])
a = np.rint(scalim)
b = a.astype(int)
im = Image.fromarray(b)
                     
im.save('/work2/jwe/stella/wifsip/20130918/science20130917A-0051EXP0009_1.png')                    
                     