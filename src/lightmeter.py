'''
Created on Apr 2, 2013

@author: jwe <jweingrill@aip.de>
'''

class LightMeter(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.time = []
        self.temp = []
        self.flux = []
        self.lux = []
        self.mag = []
        
    def loadfile(self, filename):
        """
        load the given text file
        1.4.2013;19:52:43;1,5;oC;41784;ok;

        """
        from numpy import log10,array
        from datetime import datetime
        
        f = open(filename)
        lines = f.readlines()
        f.close()
        # skip head
        lines = lines[3:]
        times = [l.split(';')[0]+' '+l.split(';')[1] for l in lines]
        self.time = [datetime.strptime(t,'%d.%m.%Y %H:%M:%S') for t in times]
        temps = [l.split(';')[2] for l in lines]
        self.temp = array([float(t.replace(',','.')) for t in temps])
        self.flux = array([float(l.split(';')[4]) for l in lines])
        self.mag = -2.5*log10(self.flux)
                
    def tolux(self, a=2.48e-7, b=0.7e-5, c=1.0e5):
        from numpy import exp, array
        flux = array(self.flux)
        return a*flux+b*exp(flux/c);
    
    def tomag(self, zero=15.2):
        from numpy import log10
        
        return -2.5*log10(self.tolux())
    
    def plot(self):
        def smooth(x):
            from numpy import ones, convolve, r_, array
            x = array(x)
            window_len=60
            s=r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
            w=ones(window_len,'d')
            y=convolve(w/w.sum(),s,mode='same')
            return y[window_len:-window_len+1]
    
        import matplotlib
        matplotlib.use('WXAgg')
        import matplotlib.pyplot as plt 
        from matplotlib.dates import HourLocator, DateFormatter
        hours = HourLocator()
        hoursFmt = DateFormatter("%H")
        #ax1 = plt.subplot(211)
        # Two subplots, the axes array is 1-d
        f, axarr = plt.subplots(2, sharex=True)
        
        axarr[0].set_ylim(7,-21)
        axarr[0].ytitle = 'flux'
        axarr[0].grid()
        axarr[0].plot_date(self.time, self.tomag(),'r')
        
        #plt.xtitle = 'time'
        axarr[1].xaxis.set_major_locator(hours)
        axarr[1].xaxis.set_major_formatter(hoursFmt)
        axarr[1].xtitle = 'time'
        axarr[1].ytitle = 'temperature (C)'
        axarr[1].grid()
        axarr[1].plot_date(self.time, smooth(self.temp),'b')
        plt.show()

if __name__ == '__main__':
    import sys
    lm = LightMeter()
    if len(sys.argv)==2:
        filename = sys.argv[1]
    else:
        filename = '/Users/jwe/Downloads/20130524T120003.txt'
    lm.loadfile(filename)
    lm.plot()