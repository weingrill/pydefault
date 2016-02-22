'''
Created on Apr 26, 2013

@author: jwe
'''
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    import astropy.io.fits as pyfits
    from functions import scaleto
    from numpy import where, clip
    
    hdulist = pyfits.open('/work2/jwe/cat/BrightStarCatalogue.fit')
        #hdr = hdulist[0].header
    tab = hdulist[1].data
    hdulist.close()
    #tab = tab[::5]
    i = where(tab['Vmag']<7.0)
    tab = tab[i]
    ra = tab['RAJ2000']
    dec = tab['DEJ2000']
    vmag = tab['Vmag']
    bv = tab['B-V']
    
    bv = clip(bv,-0.3,1.52) # B0 to M4
    
    pts = scaleto(vmag,[80.0,1.0])
    print min(vmag),max(vmag),len(vmag)
    #plt.set_cmap('Spectral')
    plt.subplot(111, axisbg='darkblue')
    #Plejades, Hyades 03 47 00 +24 07.0
    
    bmkfov = 7.25 # degrees
    xc, yc = 66.6, 16.8
    x0 = xc - bmkfov/2
    x1 = xc + bmkfov/2
    y0 = yc - bmkfov/2
    y1 = yc + bmkfov/2
    
    plt.xlim(75.,50.)
    plt.ylim(10.,30.)
    plt.xlabel('R.A. (deg)')
    plt.ylabel('Dec. (deg)')
    plt.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0], 'w', linestyle='--')
    plt.grid(color='w')
    plt.scatter(ra, dec, s=pts, c=-bv, edgecolor='none', cmap=plt.cm.Spectral)  # @UndefinedVariable
    plt.show()