'''
Created on 31.08.2013

@author: jwe
'''
def ccdlist(images, ccdtype=None, names=False, long_format=False):
    """List CCD processing information"""
    import pyfits
    for filename in images:
        hdulist = pyfits.open(filename)
        hdr = hdulist[1].header
        hdulist.close()
        
        def getkey(hdr, key, alt):
            result = ''
            try:
                result =  hdr[key]
            except KeyError:
                result = alt
            return result
        
        hdrf = getkey(hdr, 'FILENAME', filename)
        hdr1 = getkey(hdr, 'NAXIS1', -1)
        hdr2 = getkey(hdr, 'NAXIS2', -1)
        hdri = getkey(hdr, 'IMGTYPE', '')
        hdro = getkey(hdr, 'Object', '')
        try:
            hdrnaxis1 =  hdr['NAXIS1']
        except KeyError:
            hdrnaxis1 = 'none'
            
        
        
        
        print '%s[%d,%d][real][%s][]:[%s]' % \
               (hdrf, hdr1, hdr2, hdri, hdro)

def hselect(images, fields, expr):
    """extract keyword values from images satisfying a selection
    expression"""
    pass

if __name__ == '__main__':
    from glob import glob
    
    images = glob('/work2/jwe/stella/wifsip/20131023/*.fitz')
    ccdlist(images)
