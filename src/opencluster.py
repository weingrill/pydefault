'''
Created on Oct 3, 2013

@author: jwe <jweingrill@aip.de>
'''
class OpenCluster(object):
    '''
    classdocs
    '''

    def __init__(self, objectname='', uname='', ra=None, dec=None, obsmode=None):
        """
        initializes the OpenCluster Object with certain values
        obsmode can be 'rot' for rotation photometry
                    or 'cmd' for BVI photometry
        """
        self.startdate = '2013-03-21'
        self.enddate = '2018-12-31'
        self.telescope = 'STELLA-I'
        self.withfocus = True
        self.withacquire = True
        self.withguiding = True
        self.title = 'SOCS'
        self.uname = uname
        self.propid = 'cluster.survey'
        self.abstract = 'Photometric monitoring of open stellar clusters'
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
                       'YSizeCCD':4100,
                       'XBinCCD':1,
                       'YBinCCD':1}
        self.sequence = {'sequence':'FullFilters',
                        'ExposureTime':30.0,
                        'ExposureRepeat':15,
                        'ExposureIncrease':'5*1,5*1.5,5*3',
                        'FilterSequence':'5*I,5*V,5*B',
                        'offset':0.0}
        self.object = {'ObjectName':objectname,
                       'RA':ra,
                       'Dec':dec}
        self.constraints = {'MoonDistance.Min':20,
                            'SolHeight.Max':-16.0,
                            'AirmassTarget.Max':2.0,
                            'AltTarget.Min':30.0}
        self.file = ''
        self.filename = ''
        self._duration = 0
        self.fields = 1

        if obsmode == 'bv':
            self.sequence['ExposureTime'] = 30.0
            self.sequence['ExposureRepeat'] = 12
            self.sequence['ExposureIncrease'] = '3*1,3*3,3*10,3*10'
            self.sequence['FilterSequence'] = '3*V,3*B,3*V,3*B'
            self.mode['mode'] = 'BlockedPerNight'
            self.mode['maxobservations'] = 20
            self.mode['pernight'] = 2 
            self.mode['timeout'] = 3600000
            self.mode['jd0'] = 2456609.0
            self.mode['zero'] = 180.0

        if obsmode == 'bvsl20':
            self.sequence['ExposureTime'] = 20.0 # 20.0 + 200.0
            self.sequence['ExposureRepeat'] = 20
            self.sequence['ExposureIncrease'] = '5*1,5*3,5*10,5*10'
            self.sequence['FilterSequence'] = '5*V,5*B,5*V,5*B'
            self.mode['mode'] = 'BlockedPerNight'
            self.mode['maxobservations'] = 20
            self.mode['pernight'] = 2 
            self.mode['timeout'] = 3600000
            self.mode['jd0'] = 2456609.0
            self.mode['zero'] = 180.0

        if obsmode == 'bvsl30':
            self.sequence['ExposureTime'] = 30.0 # 30.0 + 300.0
            self.sequence['ExposureRepeat'] = 20
            self.sequence['ExposureIncrease'] = '5*1,5*3,5*10,5*10'
            self.sequence['FilterSequence'] = '5*V,5*B,5*V,5*B'
            self.mode['mode'] = 'BlockedPerNight'
            self.mode['maxobservations'] = 20
            self.mode['pernight'] = 2 
            self.mode['timeout'] = 3600000
            self.mode['jd0'] = 2456609.0
            self.mode['zero'] = 180.0

        if obsmode == 'bvi':
            self.sequence['ExposureTime'] = 200.0 
            self.sequence['ExposureRepeat'] = 15
            self.sequence['ExposureIncrease'] = '5*1,5*1.5,5*3'
            self.sequence['FilterSequence'] = '5*I,5*V,5*B'
            self.mode['mode'] = 'BlockedPerNight'
            self.mode['maxobservations'] = 20
            self.mode['pernight'] = 2 
            self.mode['timeout'] = 3600000
            self.mode['jd0'] = 2456609.0
            self.mode['zero'] = 180.0

        if obsmode == 'rot':
            self.sequence['ExposureTime'] = 30.0
            self.sequence['ExposureRepeat'] = 2
            self.sequence['ExposureIncrease'] = '1,10'
            self.sequence['FilterSequence'] = 'V,V'
            self.mode['mode'] = 'FewPerNight'
            self.mode['timeout'] = 1320000 # duration*fields*1000
            self.mode['pernight'] = 10 # can be refined
            self.mode['impact'] = 0.5
        
        
        
        if self.object['RA'] is None or self.object['Dec'] is None:
            self.get_coordiantes()
        
    def get_coordiantes(self):
        """
        queries the object coordinates
        """
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
        """
        returns new subframes
        """
        d = 1320.2/3600.0 # 1320.2arcsec is the fov of WiFSIP
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
        
        subframes = []
        for f in fields:
            i = fields.index(f)
            subframes.append(OpenCluster(objectname=self.object['ObjectName'],
                                         uname = '%s %s' % (self.uname, names[i]), 
                                         ra = f[0], dec = f[1]))
        #do some deep copy of the object:
        for sf in subframes:
            sf.startdate = self.startdate
            sf.enddate = self.enddate
            sf.telescope = self.telescope
            sf.withfocus = self.withfocus
            sf.withacquire = self.withacquire
            sf.withguiding = self.withguiding
            sf.title = self.title
            #sf.uname = self.uname # is set at init
            sf.propid = self.propid
            sf.abstract = self.abstract
            sf.pi = self.pi
            sf.affil = self.affil
            sf.team = self.team
            sf.title = self.title
            sf.mode = self.mode
            sf.camera = self.camera
            sf.sequence = self.sequence
            #sf.object = self.object # new coordinates and object mustn't change!
            sf.constraints = self.constraints
            
        return subframes
            
    def plot(self, axis=None):
        """
        plot the fov
        """
        import matplotlib.pyplot as plt
        ra, dec = self.object['RA'],self.object['Dec']
        d2 = 0.5*1320.2/3600.0
        ras = [ra-d2, ra+d2, ra+d2, ra-d2, ra-d2]
        das = [dec-d2, dec-d2, dec+d2, dec+d2, dec-d2]
        if axis is None:
            plt.plot(ras, das)
        else:
            axis.plot(ras, das)

    @property
    def duration(self):
        """
        getter for duration
        """
        expt = float(self.sequence['ExposureTime'])
        seq = self.sequence['ExposureIncrease'] #'5*1,5*1.5,5*3'
        sseq = seq.split(',')
        repeat = 0
        duration = 0.0
        for s in sseq:
            if s.find('*')>0:
                repeat += int(s.split('*')[0])
                duration += expt*int(s.split('*')[0])*float(s.split('*')[1])
            else:
                repeat += 1
                duration += expt*float(s)
                
        self.sequence['ExposureRepeat'] = repeat
        return duration
        
    def tofile(self, path = './'):
        """
        stores the data to a save file
        """
        #ngc6940_sw_bvi_20130926.save
        import os
        from datetime import date, datetime
        self.file = self.uname.lower().replace(' ','_')+'.xml'
        filename, ext =  os.path.splitext(self.file)
        todaystr = date.strftime(date.today(),'%Y%m%d')
        self.filename = os.path.join(path,filename + '_' + todaystr + '.save')
        def bin(boolean):
            '''return the binary value to an appropriate string'''
            return str(boolean).lower()
        
        f = open(self.filename, 'wt')
        f.write('#User input\n')
        #Thu Sep 26 10:33:48 CEST 2013
        f.write('#%s\n' % datetime.strftime(datetime.now(),'%a %b %d %H:%M:%S %Z %Y'))
        f.write('startdate=%s\n' % self.startdate)
        f.write('enddate=%s\n' % self.enddate)
        f.write('TELESCOPE=%s\n' % self.telescope)
        f.write('withfocus=%s\n' % bin(self.withfocus))
        f.write('withacquire=%s\n' % bin(self.withacquire))
        f.write('withguiding=%s\n' % bin(self.withguiding))
        f.write('title=%s\n' % str(self.title))
        f.write('uname=%s\n' % str(self.uname))
        f.write('propid=%s\n' % self.propid)
        f.write('abstract=%s\n' % self.abstract)
        f.write('pi=%s\n' % self.pi)
        f.write('affil=%s\n' % self.affil)
        f.write('team=%s\n' % self.team)
        f.write('mode=%s\n' % self.mode['mode'])
        if self.mode['mode']=='FewPerNight':
            f.write('mode.timeout=%d\n' % self.mode['timeout'])
            f.write('mode.pernight=%d\n' % self.mode['pernight'])
            f.write('mode.impact=%.1f\n' % self.mode['impact'])
            
        if self.mode['mode']=='BlockedPerNight':
            f.write('mode.maxobservations=%d\n' % self.mode['maxobservations'])
            f.write('mode.pernight=%d\n' % self.mode['pernight'])
            f.write('mode.timeout=%d\n' % self.mode['timeout'])
            f.write('mode.jd0=%.1f\n' % self.mode['jd0'])
            f.write('mode.zero=%.1f\n' % self.mode['zero'])
            
        f.write('camera=%s\n' % self.camera['camera'])
        f.write('camera.XOffCCD=%d\n' % self.camera['XOffCCD'])
        f.write('camera.YOffCCD=%d\n' % self.camera['YOffCCD'])
        f.write('camera.XSizeCCD=%d\n' % self.camera['XSizeCCD'])
        f.write('camera.YSizeCCD=%d\n' % self.camera['YSizeCCD'])
        f.write('camera.XBinCCD=%d\n' % self.camera['XBinCCD'])
        f.write('camera.YBinCCD=%d\n' % self.camera['YBinCCD'])
        f.write('sequence=%s\n' % self.sequence['sequence'])
        f.write('sequence.ExposureTime=%d\n' % self.sequence['ExposureTime'])
        f.write('sequence.ExposureRepeat=%d\n' % self.sequence['ExposureRepeat'])
        f.write('sequence.ExposureIncrease=%s\n' % self.sequence['ExposureIncrease'])
        f.write('sequence.FilterSequence=%s\n' % self.sequence['FilterSequence'])
        f.write('sequence.offset=%.1f\n' % self.sequence['offset'])
        f.write('object.ObjectName=%s\n' % self.object['ObjectName'])
        f.write('object.RA=%.5f\n' % self.object['RA'])
        f.write('object.Dec=%.6f\n' % self.object['Dec'])
        f.write('constraints.MoonDistance.Min=%d\n' % self.constraints['MoonDistance.Min'])
        f.write('constraints.SolHeight.Max=%.1f\n' % self.constraints['SolHeight.Max'])
        f.write('constraints.AirmassTarget.Max=%.1f\n' % self.constraints['AirmassTarget.Max'])
        f.write('constraints.AltTarget.Min=%.1f\n' % self.constraints['AltTarget.Min'])
        f.write('file=%s\n' % self.file)
        f.write('duration=%d\n' % self.duration)
        f.flush()
        f.close()

        fin = open(self.filename, 'rt')
        lines = fin.readlines()
        fin.close()

    def fromfile(self):
        """
        loads the data from a save file
        """
        f = open(filename, 'rt')
        lines = f.readlines()
        f.close
        #TODO: processing of lines
        
    def transfer(self, target='sro@stella:/stella/home/www/uploads/weingrill/save/'):
        """
        uploads the files to stella for the submission tool
        """
        
        from subprocess import call
        import os
        source = self.filename
        call(['/usr/bin/scp', source, target])

    def tycho(self):
        """
        plot tycho stars upon fov
        """
        from tycho import tycho
        
        tycho(self.object['RA'], 
              self.object['Dec'],
              fov = 0.2*1.4142, 
              grid=False, 
              background=False, 
              show=False)

def do_ngc2281():
    import time

    ngc2281 = OpenCluster(objectname='NGC 2281', uname='NGC 2281 rot', obsmode='rot')           
    ngc2281_subframes = ngc2281.plan_wifsip(nfields=5)
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for sf in ngc2281_subframes:
        print sf.uname
        sf.plot()
        sf.tycho()
        time.sleep(1) # otherwise submit.jnlp gets confused
        sf.tofile('/home/jwe/')
        sf.transfer()
    xmin,xmax = plt.xlim()
    plt.xlim(xmax,xmin)
    plt.show()

def do_ngc1674():
    import time

    ngc2281 = OpenCluster(objectname='NGC 1647', uname='NGC 1647 bv', obsmode='bv')           
    ngc2281_subframes = ngc2281.plan_wifsip(nfields=5)
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1.)

    for sf in ngc2281_subframes:
        print sf.uname
        sf.plot(ax)
        sf.tycho()
        time.sleep(1) # otherwise submit.jnlp gets confused
        sf.tofile('/home/jwe/')
        sf.transfer()
    xmin,xmax = plt.xlim()
    ax.set_xlim(xmax,xmin)
    plt.show()

def do_ngc2281center():
    import time

    ngc2281 = OpenCluster(objectname='NGC 2281', uname='NGC 2281 BV center', obsmode='bvsl20')           
    ngc2281.plot()
    ngc2281.tycho()
    time.sleep(1) # otherwise submit.jnlp gets confused
    ngc2281.tofile('/home/jwe/')
    ngc2281.transfer()
    xmin,xmax = plt.xlim()
    plt.xlim(xmax,xmin)
    plt.show()
        
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt
    #do_ngc2281center()
    do_ngc1674()    
    
    