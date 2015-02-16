'''
Created on 11.02.2015

@author: Joerg Weingrill <jweingrill@gmail.com>
'''

if __name__ == '__main__':
    import ephem
    
    kpno = ephem.Observer()
    kpno.date = '2015/02/11 18:00:00'
    kpno.lon, kpno.lat = '-111.600578', '31.958036'
    kpno.horizon = '20'
    kpno.elevation = 2096

    m48 = ephem.readdb("M48,f|M|F7,08:13:54.9,-05:44:06.0,5.8,2000")
    
    m48.compute(kpno)
    rise = ephem.localtime(kpno.next_rising(m48, use_center=True))
    transit = ephem.localtime(kpno.next_transit(m48))
    m48_set = ephem.localtime(kpno.next_setting(m48))

    sun = ephem.Sun()  # @UndefinedVariable
    sunset =kpno.next_setting(sun) #Sunset

    #We relocate the horizon to get twilight times
    kpno.horizon = '-18' #-6=civil twilight, -12=nautical, -18=astronomical
    beg_twilight=ephem.localtime(kpno.next_setting(sun, use_center=True)) #Begin civil twilight
    end_twilight=ephem.localtime(kpno.next_rising(sun, use_center=True)) #End civil twilight
    print 'dark     ', beg_twilight    
    print 'rise:    ', rise
    print 'transit: ', transit
    print 'set:     ', m48_set
    print 'end dark ', end_twilight    
    
