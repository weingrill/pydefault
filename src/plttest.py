'''
Created on Apr 10, 2013

@author: jwe <jweingrill@aip.de>
'''

from astronomy import airmass, height
    
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
from numpy import linspace, pi, arcsin, array

def update_ax2(ax1):
    y1, y2 = ax1.get_ylim()
    yticks = [1.0,1.1,1.3,1.5,2.0,2.5]
    ax2.set_yticks(yticks)
    ax2.set_ylim(airmass(y1), airmass(y2))
    ax2.figure.canvas.draw()


a = 3.0
print 'the height of airmass %.1f is %.2f' % (a, height(a))

fig, ax1 = plt.subplots()  # ax1 is the height scale
ax2 = ax1.twinx()          # ax2 is the airmass scale


# automatically update ylim of ax2 when ylim of ax1 changes.
ax1.callbacks.connect("ylim_changed", update_ax2)
t = linspace(0, 1.0)
y = 90.0*t**2.0
ax1.plot(t, y)

ax1.set_title('Two scales: height and airmass')
ax1.set_ylabel('height')
ax2.set_ylabel('airmass')
#airmasses = array([1.0,1.1,1.3,1.5,2.0,2.5,2.8])
#heights = [height(am) for am in airmasses]
#yticks = [height(am) for am in airmasses]
#yticklabels = [str(round(am,1)) for am in airmasses]
#ax2.set_yticks(yticks)
#ax2.set_yticklabels(yticklabels)
plt.show()

"""
if __name__ == '__main__':
    from astronomy import airmass
    
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as plt
    
    from numpy import linspace, pi, arcsin, arccos
    h = arcsin(linspace(0, 1.0))*180./pi
    a = airmass(h)
    print airmass([30., 90.])
    plt.xlabel('airmass')
    plt.ylabel('height') 
    plt.grid()
    #plt.plot(a,h)
    t = linspace(0, 1.0)
    y = arcsin(t)*180./pi
    plt.plot(t, y)
    plt.show()

"""