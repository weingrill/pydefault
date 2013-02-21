'''
Created on Feb 8, 2013

@author: jwe
'''
from astronomy import dd2hms, dd2dms

f = open('/work1/jwe/exoplanets/hoststars.tst')
lines = f.readlines()
f.close()
lines = lines[18:-1]
for l in lines:
    obj, ra, dec = l.split('\t')
    h,m,s = dd2hms(float(ra))
    dd,dm,ds = dd2dms(float(dec))
    obj = obj.replace(' ','')
    print '%-20s %02d %02d %04.1f %02d %02d %02.0f 2000' % (obj, h,m,s, dd,dm,ds)
