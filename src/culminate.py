'''
Created on Nov 14, 2013

@author: jweingrill@aip.de

Calculates the culmination of an object in Tenerife

http://www.iceinspace.com.au/forum/showthread.php?t=66580

The approx. RA at midnight on the 21st of the month is:

Month, RA
1, 8
2, 10
3, 12 [equinox]
4, 14
5, 16
6, 18
7, 20
8, 22
9, 0 [equinox]
10, 2
11, 4
12, 6

The RA on the meridian at 6 pm on the 21st of the month is equal to the month number times two.
At 8 pm RA = 2M + 2 e.g. For Jan, RA = 2x1+2 = 4 at 8 pm 
At midnt RA = 2M + 6 e.g. For Jan, RA = 2x1+6 = 8 at midnt 
At 4 am RA = 2M + 10 e.g. For Jan, RA = 2x1+10 = 12 at 4 am
'''
if __name__ == '__main__':
    #import matplotlib
    #matplotlib.use('WXAgg')
    #import matplotlib.pyplot as plt

    import ephem
    
    izana = ephem.Observer()
    date = ephem.Date('2013/12/31 12:00:00')
    izana.date = '2013/03/15'
    izana.lat = '28.301195'
    izana.lon = '-16.509209'
    sun = ephem.Sun()  # @UndefinedVariable
        
    for i in range(365):
        newdate = date + i
        izana.date = newdate
        sun.compute(izana)
        at = izana.next_antitransit(sun)
        
        izana.date = at 
        print izana.sidereal_time(), izana.date 
    
