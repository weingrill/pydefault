'''
Created on Apr 26, 2013

@author: jwe
'''
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt

    import pyfits
    import pylab
    from functions import scaleto
    from numpy import where, clip
    
    hdulist = pyfits.open('/work2/jwe/cat/BrightStarCatalogue.fit')
        #hdr = hdulist[0].header
    tab = hdulist[1].data
    hdulist.close()
    #tab = tab[::5]
    i = where(tab['Vmag']<6.5)
    tab = tab[i]
    ra = tab['RAJ2000']
    dec = tab['DEJ2000']
    vmag = tab['Vmag']
    bv = tab['B-V']
    
    bv = clip(bv,-0.3,1.52) # B0 to M4
    
    pts = scaleto(vmag,[40.0,1.0])
    print min(vmag),max(vmag),len(vmag)
    #plt.set_cmap('Spectral')
    plt.subplot(111, axisbg='darkblue')
    #Plejades, Hyades 03 47 00 +24 07.0
    plt.xlim(75.,50.)
    plt.ylim(10.,30.)
    plt.grid()
    plt.scatter(ra, dec, s=pts, c=-bv, edgecolor='none', cmap=pylab.cm.Spectral)
    plt.show()