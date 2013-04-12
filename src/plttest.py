'''
Created on Apr 10, 2013

@author: jwe <jweingrill@aip.de>
'''

if __name__ == '__main__':
    from astronomy import airmass
    
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt
    
    from numpy import arange
    h = arange(91)
    a = airmass(h)
    print airmass([30., 90.])
    plt.xlabel('airmass')
    plt.ylabel('height') 
    plt.grid()
    plt.plot(a,h)
    plt.show()
