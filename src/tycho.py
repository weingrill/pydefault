'''
Created on Oct 17, 2013

@author: jwe
'''
def tycho(ra_center, dec_center, fov = 1.0, grid=True, background=True, show=True):
    import matplotlib.pyplot as plt

    import pylab
    from functions import scaleto
    from numpy import clip, array
    from datasource import DataSource
    from pylab.cm import Spectral  # @UnresolvedImport
    
    tycho = DataSource(database='wifsip', host='pina', user='sro')
    query = """select mradeg, mdedeg, vt, bt-vt
         from tycho where circle(point(%f,%f),%f) @> circle(coord,0)
         and not bt-vt is null;""" % (ra_center, dec_center, fov)
    data = tycho.query(query)
    ra = array([float(d[0]) for d in data])
    dec = array([float(d[1]) for d in data])
    vmag = array([float(d[2]) for d in data])
    bv = array([float(d[3]) for d in data])
    tycho.close()
        
    bv = clip(bv,-0.3,1.52) # B0 to M4
    
    pts = scaleto(vmag,[80.0,1.0])
    #print min(vmag),max(vmag),len(vmag)
    if background:
        plt.subplot(111, axisbg='darkblue')
    #Plejades, Hyades 03 47 00 +24 07.0
    if grid:
        plt.grid(color='w')
    plt.scatter(ra, dec)
    plt.scatter(ra, dec, s=pts, c=-bv, edgecolor='none', cmap=Spectral)
    if show:
        plt.show()
    
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1.)
    tycho(114.145833333, -14.483, show=False) #-14.4833333333
    xmin,xmax = plt.xlim()
    plt.xlim(xmax,xmin)
    plt.show()
