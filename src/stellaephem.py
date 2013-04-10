'''
Created on Apr 10, 2013

@author: jwe <jweingrill@aip.de>
'''

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.transforms as mtransforms
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
    
    import ephem
    import datetime
    from numpy import pi,empty, sin
    
    stella = ephem.Observer()
    #stella.lon, stella.lat = '-16.50925', '28.301215'
    stella.lon, stella.lat = '13.104659', '52.404963'
    sun, moon = ephem.Sun(), ephem.Moon()
    #stella.date = '2013/4/10 12:00:00'
    stella.pressure = 0
    stella.horizon = '-0:34'
    stella.elevation = 2000
    #sun.compute(stella)
    #moon.compute(stella)
    print 'Moonrise:', stella.previous_rising(moon)
    print 'Moonset: ', stella.next_setting(moon)
    print 'Sunrise: ', stella.previous_rising(sun)
    print 'Sunset:  ', stella.next_setting(sun)

    today = datetime.datetime.today()
    dt =  datetime.timedelta(days=14)
    today += dt
    
    sun_alt = empty(24)
    moon_alt = empty(24)
    hours = range(24)
    for h in hours:
        today = today.replace(hour=h,minute=0,second=0)
        stella.date = today 
        sun.compute(stella)
        moon.compute(stella)
        sun_alt[h] = float(sun.alt)
        moon_alt[h] = float(moon.alt)
        #print alt[h]
    
    fig = plt.figure()
    ax_h = SubplotHost(fig, 1,1,1)
    
    h_to_airmass = 1.0/sin((h + 244/(165 + 47*h**1.1))*pi/180)
    
    ax_h.plot()
    
    aux_trans = mtransforms.Affine2D().scale(1., h_to_airmass)
    ax_airmass = ax_h.twin(aux_trans)
    ax_airmass.set_viewlim_mode("transform")
    
    fig.add_subplot(ax_h)
    
    ax_h.set_xticks(hours)
    ax_h.grid()
    ax_h.plot(hours, sun_alt*180.0/pi,'yo')
    ax_h.plot(hours, moon_alt*180.0/pi,'go')
    
    ax_h.axis["bottom"].set_label("hours")
    ax_h.axis["left"].set_label("height (degrees)")
    ax_airmass.axis["right"].set_label("airmass")
    #ax_airmass.axis["top"].label.set_visible(True)
    ax_airmass.axis["top"].major_ticklabels.set_visible(False)
    
    ax_h.set_ylim(0,90) 
    ax_h.set_xlim(0,24) 
    
    plt.draw()  
    plt.show()    