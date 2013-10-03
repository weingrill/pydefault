'''
Created on Oct 3, 2013

@author: jwe <jweingrill@aip.de>
'''
class OpenCluster(object):
    '''
    classdocs
    '''

    def __init__(self, name=''):
        self.name = name
        self.get_coordiantes()
    
    def get_coordiantes(self):
        import PySimbad
        import astronomy as ast
        self.coords = PySimbad.SimbadCoord(self.name)
        
        if self.coords.find('+')>0:
            p = self.coords.find('+')
        else:
            p = self.coords.find('-')
        ra = self.coords[:p].strip()    
        dec = self.coords[p:].strip()
        self.ra = ast.hms2dd(ra)
        self.dec = ast.dms2dd(dec)
        self.coords=[self.ra,self.dec]
        
    def plan_wifsip(self, nfields=4):
        d = 1320.2/3600.0
        d2 = d/2
        cra, cdec = self.coords
        if nfields == 4:
            names = ['NW','NE','SW','SE']
            fields = [(cra - d2, cdec + d2),
                 (cra + d2, cdec + d2),
                 (cra - d2, cdec - d2),
                 (cra + d2, cdec - d2)]
        if nfields == 5:
            names = ['center','NW','NE','SW','SE']
            fields = [(cra    , cdec),
                 (cra - d2, cdec + d2),
                 (cra + d2, cdec + d2),
                 (cra - d2, cdec - d2),
                 (cra + d2, cdec - d2)]
        
        for f in fields:
            i = fields.index(f)
            print "%s_%s: %.5f, %.5f" % (self.name, names[i],f[0],f[1])
            
            
    def wifsip_fov(self):
        ra = coords[0]
        dec = coords[1]
        d = 1320.2/3600.0
        ras = [ra-d/2, ra+d/2, ra+d/2, ra-d/2, ra-d/2]
        das = [dec-d/2, dec-d/2, dec+d/2, dec+d/2, dec-d/2]
        plt.plot(ras, das, color)
        
        
if __name__ == '__main__':
    ngc2281 = OpenCluster('NGC 2281') #1850
    ngc2281.get_coordiantes()
    print ngc2281.coords
    ngc2281.plan_wifsip(nfields=4)
    
    ic4756 = OpenCluster('IC 4756')
    ic4756.get_coordiantes()
    print ic4756.ra,ic4756.dec

