'''
Created on May 14, 2014

@author: jwe
'''

class TimeLine(object):
    '''
    produces a timeline of observations for a given object
    '''


    def __init__(self, target):
        '''
        Constructor
        '''
        from datasource import DataSource
        import datetime
        
        self.target = target
        
        wifsip = DataSource(database='wifsip', 
                           user='sro', 
                           host='pina.aip.de')
        
        query = """SELECT datesend, expt, 1.*matched/stars
        FROM frames 
        WHERE object like '%s %%' 
        ORDER BY datesend;""" % (self.target)
        
        result = wifsip.query(query)
        
        self.dates = [r[0] for r in result]
        self.expt = [r[1] for r in result]
        self.stars = [r[2] for r in result]
        
        self.end = self.dates
        for i in range(len(self.dates)):
            self.end[i] = self.dates[i] + \
            datetime.timedelta(seconds=self.expt[i])
        
        
    def plot(self):
        import matplotlib.pyplot as plt
        
        fig = plt.figure(figsize=(10,6))
        plt.plot_date(x=self.dates, y=self.stars)
        fig.autofmt_xdate()
        plt.title('%s' % self.target)
        plt.grid(which='both')
        plt.ylabel('fraction of matched stars')
        plt.savefig('/home/jwe/Downloads/%s timeline.pdf' % self.target)
        #plt.show()
        plt.close()
    
if __name__ == '__main__':
    tl = TimeLine("NGC 1647")
    tl.plot()