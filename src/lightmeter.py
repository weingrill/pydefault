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
        from numpy import log10
        from datetime import datetime
        
        f = open(filename)
        lines = f.readlines()
        f.close()
        # skip head
        lines = lines[3:]
        for line in lines:
            cols = line.split(';')
            if cols[0] != '#':
                t = cols[0].split('.')
                dd,mo,yy = int(t[0]), int(t[1]), int(t[2])
                
                t = cols[1].split(':')
                hh,mm,ss = int(t[0]), int(t[1]), int(t[2])
                self.time.append(datetime(yy,mo,dd,hh,mm,ss))
                self.temp.append(float(cols[2].replace(',','.')))
                flux = float(cols[4])
                self.flux.append(flux)
                
                self.mag.append(-2.5*log10(float(cols[4])))

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
        ax1 = plt.subplot(211)
        ax1.xaxis.set_major_locator(hours)
        ax1.xaxis.set_major_formatter(hoursFmt)
        plt.ylim(7,-21)
        plt.xtitle = 'time'
        plt.ytitle = 'flux'
        plt.grid()
        plt.plot_date(self.time, self.tomag(),'r')
        ax2 = plt.subplot(212)
        ax2.xaxis.set_major_locator(hours)
        ax2.xaxis.set_major_formatter(hoursFmt)
        plt.xtitle = 'time'
        plt.ytitle = 'temperature (C)'
        plt.grid()
        plt.plot_date(self.time, smooth(self.temp),'b')
        plt.show()

if __name__ == '__main__':
    lm = LightMeter()
    lm.loadfile('/home/jwe/Downloads/20130403T120002.txt')
    lm.plot()