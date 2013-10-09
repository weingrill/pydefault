'''
Created on Oct 3, 2013

@author: jwe <jweingrill@aip.de>
'''
class OpenCluster(object):
    '''
    classdocs
    '''

    def __init__(self, objectname='', uname=''):
        self.startdate = '2013-03-21'
        self.enddate = '2014-12-31'
        self.telescope = 'STELLA-I'
        self.withfocus = True
        self.withacquire = True
        self.withguiding = True
        self.title = ''
        self.uname = uname
        self.propid = 'cluster.survey'
        self.abstract = ''
        self.pi = 'Weingrill'
        self.affil = 'CORE'
        self.team = 'Barnes'
        self.mode = {'mode':'FewPerNight',
                     'timeout':3600000,
                     'pernight':2,
                     'impact':1.0}
        self.camera = {'camera':'direct',
                       'XOffCCD':0,
                       'YOffCCD':0,
                       'XSizeCCD':2050,
                       'YSizeCCD':2050,
                       'XBinCCD':1,
                       'YBinCCD':1}
        self.sequence = {'sequence':'FullFilters',
                        'ExposureTime':200.0,
                        'ExposureRepeat':15,
                        'ExposureIncrease':'5*1,5*1.5,5*3',
                        'FilterSequence':'5*I,5*V,5*B',
                        'offset':0.0}
        self.object = {'ObjectName':objectname,
                       'RA':0.0,
                       'Dec':0.0}
        self.constraints = {'MoonDistance.Min':45,
                            'SolHeight.Max':-16.0,
                            'AirmassTarget.Max':2.0,
                            'AltTarget.Min':30.0}
        self.file = ''
        self.duration = 0
        self.get_coordiantes()
        
    def get_coordiantes(self):
        import PySimbad
        import astronomy as ast
        
        self.coords = PySimbad.SimbadCoord(self.object['ObjectName'])
        
        if self.coords.find('+')>0:
            p = self.coords.find('+')
        else:
            p = self.coords.find('-')
        ra = self.coords[:p].strip()    
        dec = self.coords[p:].strip()
        self.object['RA'] = ast.hms2dd(ra)
        self.object['Dec'] = ast.dms2dd(dec)
        self.coords=[self.object['RA'],self.object['Dec']]
        
    def plan_wifsip(self, nfields=4):
        d = 1320.2/3600.0
        d2 = d/2
        cra, cdec = self.object['RA'],self.object['Dec']
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
            print "%s_%s: %.5f, %.5f" % (self.object['ObjectName'], names[i],f[0],f[1])
            
            
    def wifsip_fov(self):
        ra, dec = self.object['RA'],self.object['Dec']
        
        d = 1320.2/3600.0
        ras = [ra-d/2, ra+d/2, ra+d/2, ra-d/2, ra-d/2]
        das = [dec-d/2, dec-d/2, dec+d/2, dec+d/2, dec-d/2]
        plt.plot(ras, das, color)
    
    def tofile(self, path = './'):
        #ngc6940_sw_bvi_20130926.save
        import os
        from datetime import date, datetime
        self.file = self.uname.lower().replace(' ','_')
        filename, ext =  os.path.splitext(self.file)
        todaystr = date.strftime(date.today(),'%Y%m%d')
        filename = os.path.join(path,filename + '_' + todaystr + '.save')
        f = open(filename, 'wt')
        f.write('#User input\n')
        #Thu Sep 26 10:33:48 CEST 2013
        nowstr = datetime.strftime(datetime.now(),'%a %b %d %H:%M:%S %Z %Y')
        f.write('#'+nowstr+'\n')
        f.write('startdate='+self.startdate+'\n')
        f.write('enddate='+self.enddate+'\n')
        f.write('TELESCOPE='+self.telescope+'\n')
        f.write('withfocus='+str(self.withfocus).lower()+'\n')
        f.write('withacquire='+str(self.withacquire).lower()+'\n')
        f.write('withguiding='+str(self.withguiding).lower()+'\n')
        f.write('title='+str(self.title)+'\n')
        f.write('uname='+str(self.uname)+'\n')
        f.write('propid='+self.propid+'\n')
        f.write('abstract='+self.abstract+'\n')
        f.write('pi='+self.pi+'\n')
        f.write('affil='+self.affil+'\n')
        f.write('team='+self.team+'\n')
        f.write('mode='+self.mode['mode']+'\n')
        f.write('mode.timeout=%d\n' % self.mode['timeout'])
        f.write('mode.pernight=%d\n' % self.mode['pernight'])
        f.write('mode.impact=%.1f\n' % self.mode['impact'])
        f.write('camera='+self.camera['camera']+'\n')
        f.write('camera.XOffCCD=%d\n' % self.camera['XOffCCD'])
        f.write('camera.YOffCCD=%d\n' % self.camera['YOffCCD'])
        f.write('camera.XSizeCCD=%d\n' % self.camera['XSizeCCD'])
        f.write('camera.YSizeCCD=%d\n' % self.camera['YSizeCCD'])
        f.write('camera.XBinCCD=%d\n' % self.camera['XBinCCD'])
        f.write('camera.YBinCCD=%d\n' % self.camera['YBinCCD'])
        f.write('sequence='+self.sequence['sequence']+'\n')
        f.write('sequence.ExposureTime=%.1f\n' % self.sequence['ExposureTime'])
        f.write('sequence.ExposureRepeat=%d\n' % self.sequence['ExposureRepeat'])
        f.write('sequence.ExposureIncrease=%s\n' % self.sequence['ExposureIncrease'])
        f.write('sequence.FilterSequence=%s\n' % self.sequence['FilterSequence'])
        f.write('sequence.offset=%.1f\n' % self.sequence['offset'])
        f.write('object.ObjectName=%s\n' % self.object['ObjectName'])
        f.write('object.RA=%.5f\n' % self.object['RA'])
        f.write('object.Dec=%.6f\n' % self.object['Dec'])
        f.write('constraints.MoonDistance.Min=%d\n' % self.constraints['MoonDistance.Min'])
        f.write('constraints.SolHeight.Max=%d\n' % self.constraints['SolHeight.Max'])
        f.write('constraints.AirmassTarget.Max=%d\n' % self.constraints['AirmassTarget.Max'])
        f.write('constraints.AltTarget.Min=%d\n' % self.constraints['AltTarget.Min'])
        f.write('file=%s\n' % self.file)
        f.write('duration=%d\n' % self.duration)
        f.close()
    
    def fromfile(self):
        pass    
        
if __name__ == '__main__':
    ngc2281 = OpenCluster(objectname='NGC 2281', uname='NGC2281 SW BVI') #1850
    ngc2281.get_coordiantes()
    ngc2281.plan_wifsip(nfields=4)
    ngc2281.tofile('/home/jwe/')
    
